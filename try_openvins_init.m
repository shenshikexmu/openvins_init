clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');
addpath('repropagate');
addpath('quaternion');

global matlab_or_octave
matlab_or_octave=0; 

    
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


global min_px_dist  grid_x grid_y num_features threshold currid

gravity_mag= 9.81;

min_px_dist=10;
grid_x=5; 
grid_y=5;
currid=0;


pyr_levels=5;

num_pts=200;  % 定位需要提取的特征点数
num_cameras=1;

num_features=round(num_pts/num_cameras);

fast_threshold=20;

threshold=fast_threshold;

win_size = [15, 15];

camK=[cam0Para.intrinsics(1),0,cam0Para.intrinsics(3);...
      0,cam0Para.intrinsics(2),cam0Para.intrinsics(4);...
      0,0,1];

camD=[cam0Para.distortion_coefficients(1),cam0Para.distortion_coefficients(2),cam0Para.distortion_coefficients(3),cam0Para.distortion_coefficients(4)];

camR=cam0Para.T_BS.data(1:3,1:3);

camT=cam0Para.T_BS.data(1:3,4);


max_focallength=max(camK(1,1),camK(2,2));


global ACC_N GYR_N ACC_W GYR_W 


ACC_N=0.5;                     %bao
GYR_N=0.1;
ACC_W=0.001;
GYR_W=0.00001;




i=989;



m=16;
for n=1:m
    
    img{n}=imread([file_cam0,'data/',datacsv_cam0{i+n-1,2}]);
    img{n}=cv_equalizeHist(img{n});
    imgpyr{n}=cv_buildOpticalFlowPyramid(img{n},win_size,pyr_levels);
    
    if (exist(['pts_ids_',num2str(i),'_',num2str(m),'.mat'])==2)   
        if n==m           
            load(['pts_ids_',num2str(i),'_',num2str(m),'.mat']);
        end 
        continue;   
    end
      
    if n==1
        mask=getMask(img{n});
        pts{n}=[];
        ids{n}=[];
        [pts{n},ids{n}]=perform_detection_monocular(imgpyr{n},mask,pts{n},ids{n});
        continue;   
    end
    
    criteria.max_iters=30;
    criteria.epsilon=0.01;
    cv_OPTFLOW_USE_INITIAL_FLOW=1;
    
    pts{n}=pts{n-1};
    [pts{n}, status, err] = cv_calcOpticalFlowPyrLK(imgpyr{n-1}, imgpyr{n}, pts{n-1}, pts{n}, win_size, pyr_levels, criteria, cv_OPTFLOW_USE_INITIAL_FLOW);

    pts{n-1}=refine(pts{n-1},status);
    ids{n-1}=refine(ids{n-1},status);
    pts{n}=refine(pts{n},status);

    pts_n{n-1}=zeros(size(pts{n-1},1),2);
    pts_n{n}=zeros(size(pts{n},1),2);

    for p=1:size(pts{n-1},1)
        pts_n{n-1}(p,:)=undistort_cv(pts{n-1}(p,1:2)-[1,1], camK,camD);
        pts_n{n}(p,:)=undistort_cv(pts{n}(p,1:2)-[1,1], camK,camD);
    end

    [mask]=cv_findFundamentalMat(pts_n{n-1}, pts_n{n}, 'cv_FM_RANSAC', 1/max_focallength ,0.999);

    pts{n-1}=refine(pts{n-1},mask);
    pts_n{n-1}=refine(pts_n{n-1},mask);
    ids{n-1}=refine(ids{n-1},mask);
    pts{n}=refine(pts{n},mask);
    pts_n{n}=refine(pts_n{n},mask);
    ids{n}=ids{n-1};
    
    if n==m
        save(['pts_ids_',num2str(i),'_',num2str(n),'.mat'],'pts','ids','pts_n');
    end
    
end

%%

ba=[0;0;0];
bg=[0;0;0];

for n=1:m

    camTimeStamp(n,1)=datacsv_cam0{i+n-1,1}*10e-10;

    if n>1

        imuData_fragment{n}=get_Imu_Fragment_tc1_tc2(imuData,camTimeStamp(n-1),camTimeStamp(n));

    else

        imuData_fragment{n}=[];

    end

    imuPropagate{n}=repropagate_VINS_Mono(imuData_fragment{n},ba,bg);

    imuData_fragment_joint{n}=get_Imu_Fragment_tc1_tc2(imuData,camTimeStamp(1),camTimeStamp(n));

    imuPropagate_joint{n}=repropagate_VINS_Mono(imuData_fragment_joint{n},ba,bg);

    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  paper <OpenVINS State Initiaolization: Details and Derivations> 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_frame=12;
size_ids=size(ids{size_frame},1);
M0=zeros(size_ids*size_frame*2,size_ids*3+6);
b0=zeros(size_ids*size_frame*2,1);

for i=1:size_ids%size(ids{10},1)

    for n=1:size_frame

        for s=1:size(ids{n},1)

             if ids{n}(s)==ids{size_frame}(i)

                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 pts_opti{n}(i,:)=pts{n}(s,:);
                 ids_opti{n}(i,:)=ids{n}(s,:);
                 pts_n_opti{n}(i,:)=pts_n{n}(s,:);
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 Gamma=[1,0,-pts_n{n}(s,1);...
                        0,1,-pts_n{n}(s,2)];

                 C_I_R=camR';

                 I0_alpha_k=imuPropagate_joint{n}{1};

                 I0_q_k=imuPropagate_joint{n}{2};

                 Delta_T=imuPropagate_joint{n}{6};

                 Ik_I0_R=quatern2rotMat(I0_q_k)';

                 % 论文中公式（28）（29）（30）存在错误，这里计算的时候做了修改
                 M((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,(i-1)*3+1:(i-1)*3+3)=Gamma*C_I_R*Ik_I0_R;

                 M((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,size_ids*3+1:size_ids*3+3)=-Gamma*C_I_R*Ik_I0_R*Delta_T;

                 M((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,size_ids*3+4:size_ids*3+6)=0.5*Gamma*C_I_R*Ik_I0_R*Delta_T*Delta_T;

                 b((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,1)=Gamma*camR'*camT+ Gamma*C_I_R*Ik_I0_R*I0_alpha_k;

             end

        end

    end

end

%%
A1=M(:,1:size_ids*3+3);
A2=M(:,size_ids*3+4:size_ids*3+6);

inv_A1A1=inv(A1'*A1);

tempMatrix=eye(size_ids*size_frame*2)-A1*inv_A1A1*A1';
D=A2'*tempMatrix*A2;

d=A2'*tempMatrix*b;

lambda=getLambda(D,d,gravity_mag);

I0_g=inv(D-lambda*eye(3))*d;

x1=-inv_A1A1*A1'*A2*I0_g+inv_A1A1*A1'*b;


x=M\b;

%%


I0_v_I0=x1(size_ids*3+1:size_ids*3+3,1);

G_I0_q=getInitQuaternionfromAcc_yfront(I0_g);
G_I0_R=quatern2rotMat(G_I0_q);

for i=1:size_ids

    I0_p_f(:,i)=x1(i*3-2:i*3);

    G_p_f(:,i)=G_I0_R*I0_p_f(:,i);

end



for n=1:size_frame

    I0_Ik_q=imuPropagate_joint{n}{2};
    I0_Ik_R=quatern2rotMat(I0_Ik_q);
    delta_T=imuPropagate_joint{n}{6};
    I0_alpha_k=imuPropagate_joint{n}{1};
    I0_bata_k=imuPropagate_joint{n}{3};

    I0_p_Ik=I0_v_I0*delta_T-0.5*I0_g*delta_T^2+I0_alpha_k;
    I0_v_Ik=I0_v_I0-I0_g*delta_T+I0_bata_k;

    G_Ik_q=quaternProd(G_I0_q,I0_Ik_q);  %   G_Ik_R=quatern2rotMat(G_Ik_q)    G_Ik_R_=G_I0_R*I0_Ik_R
    G_p_Ik=G_I0_R*I0_p_Ik;
    G_v_Ik=G_I0_R*I0_v_Ik;


     x_I_k(:,n)=[G_Ik_q;G_p_Ik;G_v_Ik;bg;ba];


end

drawProjection(x_I_k,G_p_f,camR,camT);


%%


[x_I_k_opti,G_p_f_opti]=optimaization_openvins_init(x_I_k,G_p_f,imuPropagate,pts_opti,pts_n_opti,camK,camD,camR,camT,gravity_mag);


%%

G_I0_q_opti=x_I_k_opti(1:4,1);
G_I0_v=x_I_k_opti(8:10,1);

I0_G_q_opti=invQuaternion(G_I0_q_opti);
I0_G_R_opti=quatern2rotMat(I0_G_q_opti);

I0_g_opti=I0_G_R_opti*[0;0;1]*9.8;

I0_v_I0=I0_G_R_opti*G_I0_v;


for i=1:size(x_I_k_opti,2)

    x_I_k_opti_temp=x_I_k_opti(:,i);

    G_Ik_q_opti=x_I_k_opti_temp(1:4);
    G_Ik_p_opti=x_I_k_opti_temp(5:7);
    G_Ik_v_opti=x_I_k_opti_temp(8:10);

    I0_Ik_q_opti=quaternProd(I0_G_q_opti,G_Ik_q_opti);
    I0_Ik_p_opti=I0_G_R_opti*G_Ik_p_opti;
    I0_Ik_v_opti=I0_G_R_opti*G_Ik_v_opti;

    I0_x_I_k_opti(:,i)=[I0_Ik_q_opti;I0_Ik_p_opti;I0_Ik_v_opti;x_I_k_opti_temp(11:16)];

end

I0_p_f_opti=I0_G_R_opti*G_p_f_opti;

x_opti=[];

for i=1:size(I0_p_f_opti,2)

    x_opti(i*3-2:i*3,1)=I0_p_f_opti(:,i);

end

x_opti=[x_opti;I0_x_I_k_opti(8:10,1);I0_g_opti];

SSS=M*x_opti-b;

SS=M*x-b;
%%

nnn=6;

P1=I0_p_f_opti(:,nnn)


vvv=camR'*(P1-camT)

vvv_norm=vvv/vvv(3);

uv_dist = distort_cv(vvv_norm(1:2), camK,camD)



pts_opti{1}(nnn,:)
                

                 

Gamma=[1,0,-pts_n_opti{1}(nnn,1);...
       0,1,-pts_n_opti{1}(nnn,2)]

Gamma*vvv


%%



for n=1:1%size_frame

    for i=1:size_ids%size(ids{10},1)

       

%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          pts_opti{n}(i,:)=pts{n}(s,:);
%          ids_opti{n}(i,:)=ids{n}(s,:);
%          pts_n_opti{n}(i,:)=pts_n{n}(s,:);
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         Gamma=[1,0,-pts_n_opti{n}(i,1);...
                0,1,-pts_n_opti{n}(i,2)];

         C_I_R=camR';

         I0_alpha_k=imuPropagate_joint{n}{1};

         I0_q_k=imuPropagate_joint{n}{2};

         Delta_T=imuPropagate_joint{n}{6};

         Ik_I0_R=quatern2rotMat(I0_q_k)';

         % 论文中公式（28）（29）（30）存在错误，这里计算的时候做了修改
         MM((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,(i-1)*3+1:(i-1)*3+3)=Gamma*C_I_R*Ik_I0_R;

         MM((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,size_ids*3+1:size_ids*3+3)=-Gamma*C_I_R*Ik_I0_R*Delta_T;

         MM((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,size_ids*3+4:size_ids*3+6)=0.5*Gamma*C_I_R*Ik_I0_R*Delta_T*Delta_T;

         bb((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,1)=Gamma*camR'*camT+ Gamma*C_I_R*Ik_I0_R*I0_alpha_k;

   

    end

end

sss=MM*x_opti-bb;


a=100;




%%

drawProjection(x_I_k,G_p_f,camR,camT);

drawProjection(x_I_k_opti,G_p_f_opti,camR,camT);

drawProjection(I0_x_I_k_opti,I0_p_f_opti,camR,camT);




%%

if 1

    for n=size_frame

        drawOpticalFlowLK(imgpyr{1}{1},imgpyr{n}{1},pts_opti{1},pts_opti{n},1,n);
          
    end

end


b=1000;







