clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');
addpath('repropagate');
addpath('quaternion');
addpath('orb_slam');



global matlab_or_octave
matlab_or_octave=1; 


file_cam0='../bag/V1_02_medium/mav0/cam0/';
file_imu0='../bag/V1_02_medium/mav0/imu0/';

if  (matlab_or_octave ==1)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  matlab   
    datacsv_cam0=readcell([file_cam0,'data.csv']);
    datacsv_imu0=readcell([file_imu0,'data.csv']);
    imuData=cell2mat(datacsv_imu0(2:end,:));
    imuData(:,1)=imuData(:,1)*10e-10;
    cam0Para = ReadYaml_matlab([file_cam0,'sensor.yaml']);            % cell2mat
    cam0Para.intrinsics=cell2mat(cam0Para.intrinsics);
    cam0Para.distortion_coefficients=cell2mat(cam0Para.distortion_coefficients);
    cam0Para.T_BS.data=reshape(cell2mat(cam0Para.T_BS.data),4,4)';
else                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  octave
    pkg load io
    pkg load image
    datacsv_cam0 = csv2cell([file_cam0,'data.csv']);
    cam0Para = readYaml_octave([file_cam0,'sensor.yaml']);
    cam0Para.T_BS.data=reshape(cam0Para.data,4,4)';
    datacsv_imu0 = csv2cell([file_imu0,'data.csv']);
    imuData=cell2mat(datacsv_imu0(2:end,:));
    imuData(:,1)=imuData(:,1)*10e-10;
end

camK=[cam0Para.intrinsics(1),0,cam0Para.intrinsics(3);...
      0,cam0Para.intrinsics(2),cam0Para.intrinsics(4);...
      0,0,1];

camD=[cam0Para.distortion_coefficients(1),cam0Para.distortion_coefficients(2),cam0Para.distortion_coefficients(3),cam0Para.distortion_coefficients(4)];

camR=cam0Para.T_BS.data(1:3,1:3);

camT=cam0Para.T_BS.data(1:3,4);

gravity_mag= 9.81;

win_size = [15, 15];

pyr_levels=5;

global ACC_N GYR_N ACC_W GYR_W 

ACC_N=2.0000e-3;                     
GYR_N=1.6968e-04;
ACC_W=3.0000e-3;
GYR_W=1.9393e-05;

features=features_get_from_txt('./data/features.txt');

map_camera_times = [];
count_valid_features=0;
num_measurements=0;
cam_id=1;

allKeys = keys(features);

for i = 1:length(allKeys)
    key = allKeys{i}; % 获取当前键

    feat = features(key);
    
    map_camera_ids(i,1)=feat.featid;
    
    if ~isempty(feat.timestamps{cam_id})
        if ~ismember(feat.timestamps{cam_id},map_camera_times)
            map_camera_times = [map_camera_times; feat.timestamps{cam_id}];
        end
    end
    
    for cam_id = 1:length(feat.uvs_norm)
        if ~isempty(feat.uvs_norm{cam_id})          
            num_measurements=num_measurements+size(feat.uvs_norm{cam_id},1);
        end
    end
    count_valid_features=count_valid_features+1;
end

map_camera_times = unique(map_camera_times);
map_camera_times=[map_camera_times,map_camera_times*0,map_camera_times*0];



 for i=1:size(map_camera_times,1)

    [n_frame,timestamps_new, min_gap]=find_image_frame_corresponding_timestamps(datacsv_cam0,map_camera_times(i,1));
    map_camera_times(i,:)=[timestamps_new,n_frame, min_gap];

end

ba=[0;0;0];
bg=[0;0;0];

for n=1:size(map_camera_times,1)

    if n>1
        imuData_fragment{n}=get_Imu_Fragment_tc1_tc2(imuData,map_camera_times(n-1,1),map_camera_times(n,1));
    else
        imuData_fragment{n}=[];
    end

    imuPropagate{n}=repropagate_VINS_Mono(imuData_fragment{n},ba,bg);

    imuData_fragment_joint{n}=get_Imu_Fragment_tc1_tc2(imuData,map_camera_times(1,1),map_camera_times(n,1));

    imuPropagate_joint{n}=repropagate_VINS_Mono(imuData_fragment_joint{n},ba,bg);


end


for n=1:size(map_camera_times,1)

    img{n}=imread([file_cam0,'data/',datacsv_cam0{map_camera_times(n,2),2}]);
    img{n}=cv_equalizeHist(img{n});
    imgpyr{n}=cv_buildOpticalFlowPyramid(img{n},win_size,pyr_levels);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  paper <Vins-mono: A robust and versatile monocular visual-inertial state estimator> 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

frame1=1;
frame2=2;

[n,pts1_n,pts2_n,pts1,pts2]=features_intersection_in_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2);

R1=eye(3);
T1=[0;0;0];

[mask,R2,T2]=cv_findFundamentalMat(pts1_n,pts2_n, 'cv_FM_RANSAC', 4/max([camK(1,1),camK(2,2)]) ,0.99999999 );



mMaxIterations=200;
mSigma=1.0;
mSigma2 = mSigma*mSigma;
mvMatches12=(1:size(pts1_n,1))';
mvSets=zeros(8,mMaxIterations);

for i=1:mMaxIterations
    mvSets(:,i)=randperm(size(pts1_n,1),8)';
end

[vbMatchesInliersF, SF, F] = FindFundamental(mvMatches12, pts1_n, pts2_n, mvSets, mMaxIterations, mSigma);

[vbMatchesInliersH, SH, H] = FindHomography(mvMatches12, pts1_n, pts2_n, mvSets, mMaxIterations, mSigma);

RH = SH/(SH+SF);

mK=eye(3);

R21=[];
t21=[];
vP3D=[];
vbTriangulated=[];

[success,R21, t21,vP3D, vbTriangulated]=ReconstructH(vbMatchesInliersH,H,mK,mvMatches12,pts1_n, pts2_n,1.0,50,mSigma2);

% if RH>0.40
%     ReconstructH(vbMatchesInliersH,H,mK,R21,t21,vP3D,vbTriangulated,1.0,50);
% else
%     ReconstructF(vbMatchesInliersF,F,mK,R21,t21,vP3D,vbTriangulated,1.0,50);
% end







%[R2,T2]=Initial_R_T(pts1_n,pts2_n);

% R2=[0.996258210477046, 0.001305172626110589, 0.08641686514500559;...
%  -0.004405289686119615, 0.9993530921139655, 0.03569300638833356;...
%  -0.08631437585857517, -0.03594014199571702, 0.9956194728476701];
% 
% T2=[0.7434306827423962;...
%  0.6672805650253606;...
%  0.04524895022661077];






features=features_p_FinA_from_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2,R1,T1,R2,T2);

%drawOpticalFlowLK_featrues(imgpyr,features,map_camera_times,cam_id,cam_id,frame1,frame2);

draw_init(features,map_camera_times,R1,T1,R2,T2,cam_id,cam_id,frame1,frame2);

x_I_k(1:4,1)=rotMat2qRichard(R1);
x_I_k(5:7,1)=T1;
x_I_k(1:4,2)=rotMat2qRichard(R2);
x_I_k(5:7,2)=T2;

for i=3:size(map_camera_times,1)

    frame=i;
    [n,pts_n,pts,wP]=features_uvNorm_wP_in_frame(features,map_camera_times,cam_id,frame);
    
    [R_temp,T_temp]=solvePnP_from_P3P_norm(wP, pts_n);
    hold on
    draw_init_camera(R_temp,T_temp,frame);

    x_I_k(1:4,i)=rotMat2qRichard(R_temp);
    x_I_k(5:7,i)=T_temp;

end
%%

[x_I_k_opti,G_p_f_opti]=optimaization_vins_mono_init_using_features(x_I_k,features,camK,camD,map_camera_times);

drawProjection(x_I_k_opti,G_p_f_opti,eye(3),zeros(3,1));

%%    

AA=[];
bb=[];

for n=3:size(map_camera_times,1)

    if n>1

        gamma_k_k_1=imuPropagate{n}{2};

        J_gamma_bw=imuPropagate{n}{7}(4:6,13:15);

        q_c_bk=x_I_k_opti(1:4,n-1);

        q_c_bk_1=x_I_k_opti(1:4,n);

        q_temp=quaternProd(quaternProd(invQuaternion(gamma_k_k_1),invQuaternion(q_c_bk)),q_c_bk_1);

        AA=[AA;0.5*J_gamma_bw];

        bb=[bb;q_temp(2:4)];

    end

end

% bg=AA\bb ;       % The calculation here is not good !!!!!! not good !!!!!



%%

n=size(map_camera_times,1);


M=zeros((n-1)*6,n*3+3+1);

b=zeros((n-1)*6,1);


R_c_b=camR';
T_b_c=camT;



for i=2:n

    alpha_k_k_1=imuPropagate{i}{1};
    beta_k_k_1=imuPropagate{i}{3};
    delta_t=imuPropagate{i}{6};

    P_c0_ck=x_I_k_opti(5:7,i-1);
    p_c0_ck_1=x_I_k_opti(5:7,i);

    R_c0_ck=quatern2rotMat(x_I_k_opti(1:4,i-1));
    R_c0_ck_1=quatern2rotMat(x_I_k_opti(1:4,i));


    M((i-2)*6+1:(i-2)*6+3,(i-2)*3+1:(i-2)*3+3)=-eye(3)*delta_t;

    M((i-2)*6+1:(i-2)*6+3,n*3+1:n*3+1)=p_c0_ck_1-P_c0_ck;

    M((i-2)*6+1:(i-2)*6+3,n*3+2:n*3+4)=0.5*eye(3)*delta_t^2;


    M((i-2)*6+4:(i-2)*6+6,(i-2)*3+1:(i-2)*3+3)=-eye(3);
   
    M((i-2)*6+4:(i-2)*6+6,(i-2)*3+4:(i-2)*3+6)= eye(3);

    M((i-2)*6+4:(i-2)*6+6,n*3+2:n*3+4)=eye(3)*delta_t;

    b((i-2)*6+1:(i-2)*6+3,1)=R_c0_ck*R_c_b*alpha_k_k_1+R_c0_ck_1*R_c_b*T_b_c-R_c0_ck*R_c_b*T_b_c;

    b((i-2)*6+4:(i-2)*6+6,1)=R_c0_ck*R_c_b*beta_k_k_1;

end

%%


A1=M(:,1:end-3);

A2=M(:,end-2:end);

inv_A1A1=inv(A1'*A1);

tempMatrix=eye((n-1)*6)-A1*inv_A1A1*A1';
D=A2'*tempMatrix*A2;

d=A2'*tempMatrix*b;

lambda=getLambda(D,d,gravity_mag);

g_c0=inv(D-lambda*eye(3))*d;

x1=-inv_A1A1*A1'*A2*g_c0+inv_A1A1*A1'*b;

scale=x1(end);

%x=M\b

q_w_c0=getInitQuaternionfromAcc_xfront(g_c0);
R_w_c0=quatern2rotMat(q_w_c0);



for i=1:size(x_I_k_opti,2)


    q_c0_ck=x_I_k_opti(1:4,i);

    P_c0_ck=x_I_k_opti(5:7,i)*scale;

    R_c0_ck=quatern2rotMat(q_c0_ck);

    
    R_w_ck=R_w_c0*R_c0_ck;

    P_w_ck=R_w_c0*P_c0_ck;

    R_w_bk=R_w_ck*R_c_b;

    P_w_bk=P_w_ck-R_w_bk*T_b_c;

    x_I_K_body(1:4,i)=rotMat2qRichard(R_w_bk);

    x_I_K_body(5:7,i)=P_w_bk;


end


G_p_f_body=R_w_c0*G_p_f_opti*scale;




drawProjection(x_I_K_body,G_p_f_body,camR,camT);











































