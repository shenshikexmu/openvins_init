clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');
addpath('repropagate');
addpath('quaternion');

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


global min_px_dist  grid_x grid_y num_features threshold currid camK camD win_size pyr_levels

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


ACC_N=2.0000e-3;                     %bao
GYR_N=1.6968e-04;
ACC_W=3.0000e-3;
GYR_W=1.9393e-05;


features=containers.Map();
cam_id=1;

i=312;

m=5;
%table_img_timestamp=zeros(m,2);

for n=1:m
    
    img{n}=imread([file_cam0,'data/',datacsv_cam0{i+n-1,2}]);
    img{n}=cv_equalizeHist(img{n});
    imgpyr{n}=cv_buildOpticalFlowPyramid(img{n},win_size,pyr_levels);
    timestamp=datacsv_cam0{i+n-1,1}*10e-10;
    %table_img_timestamp(n,:)=[n,timestamp];
    mask{n}=getMask(img{n});
    
    if (exist(['features_',num2str(i),'_',num2str(m),'.mat'])==2)   
        if n==m           
            load(['features_',num2str(i),'_',num2str(m),'.mat']);
        end 
        continue;   
    end
      
    if n==1
        
        pts{n}=[];
        ids{n}=[];
        [pts{n},ids{n}]=perform_detection_monocular(imgpyr{n},mask{n},pts{n},ids{n});
  
    else
        
        [pts{n-1},ids{n-1}]=perform_detection_monocular(imgpyr{n-1}, mask{n-1}, pts{n-1},ids{n-1});

        pts{n}=pts{n-1};
        ids{n}=ids{n-1};

        [pts{n-1}, pts{n}, mask_out]=perform_matching(imgpyr{n-1}, imgpyr{n},  pts{n-1}, pts{n}, 0, 0);

        pts{n}=refine(pts{n},mask_out);
        ids{n}=refine(ids{n},mask_out);
        pts_last_plot=refine(pts{n-1},mask_out);

        %drawOpticalFlowLK(imgpyr{n-1}{1},imgpyr{n}{1},pts_last_plot,pts{n},n-1,n);
        a=10;

    end

    for s=1:size(pts{n},1)
        npt_l=undistort_cv(pts{n}(s,1:2)-[1,1], camK,camD);   
        features=update_feature(features, ids{n}(s,1), timestamp, cam_id, pts{n}(s,1), pts{n}(s,2), npt_l(1,1), npt_l(1,2));
    end

    if n==m
        save(['features_',num2str(i),'_',num2str(n),'.mat'],'features');
    end

end

%%

features=features_eliminate_1point(features);

% frame1=9;
% frame2=16;
% 
% drawOpticalFlowLK_featrues(imgpyr,features,table_img_timestamp,cam_id,cam_id,frame1,frame2);

%%

map_camera_times = [];
count_valid_features=0;
num_measurements=0;
cam_id=1;

all_ids = keys(features);

for i = 1:length(all_ids)
    id = all_ids{i}; % 获取当前键

    feat = features(id);
    
    if ~isempty(feat.timestamps{cam_id})
        for s=1:size(feat.timestamps{cam_id},1)
            if ~ismember(feat.timestamps{cam_id}(s,1),map_camera_times)
                map_camera_times = [map_camera_times; feat.timestamps{cam_id}(s,1)];
            end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  paper <OpenVINS State Initiaolization: Details and Derivations> 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_frames=size(map_camera_times,1);
num_features=count_valid_features;
M=zeros(num_measurements*2,num_features*3+6);
b=zeros(num_measurements*2,1);

num_measurements=0;

all_ids = keys(features);

for i = 1:length(all_ids)
    id = all_ids{i}; % 获取当前键

    feat = features(id);
    
    featid=feat.featid;

    map_features_ids(i,1)=featid;
    
    for cam_id = 1:length(feat.uvs_norm)
        uvs_norm=feat.uvs_norm{cam_id};
        uvs=feat.uvs{cam_id};
        
        timestamps=feat.timestamps{cam_id};
        
        for j=1:size(timestamps,1)
            
            num_measurements=num_measurements+1;
     
            camtime=timestamps(j,1);
            
            n=find_cam_n_from_map_camera_times(camtime,map_camera_times);
            
            uv_norm=uvs_norm(j,:);
            uv=uvs(j,:);
            
            Gamma=[1,0,-uv_norm(1);...
                0,1,-uv_norm(2)];

            C_I_R=camR';

            I0_alpha_k=imuPropagate_joint{n}{1};

            I0_q_k=imuPropagate_joint{n}{2};

            Delta_T=imuPropagate_joint{n}{6};

            Ik_I0_R=quatern2rotMat(I0_q_k)';

            M(num_measurements*2-1:num_measurements*2,i*3-2:i*3)=Gamma*C_I_R*Ik_I0_R;

            M(num_measurements*2-1:num_measurements*2,num_features*3+1:num_features*3+3)=-Gamma*C_I_R*Ik_I0_R*Delta_T;

            M(num_measurements*2-1:num_measurements*2,num_features*3+4:num_features*3+6)=0.5*Gamma*C_I_R*Ik_I0_R*Delta_T*Delta_T;

            b(num_measurements*2-1:num_measurements*2,1)=Gamma*camR'*camT+ Gamma*C_I_R*Ik_I0_R*I0_alpha_k;
            
        end 
        
    end

end

a=10;

A1=M(:,1:num_features*3+3);
A2=M(:,num_features*3+4:num_features*3+6);

inv_A1A1=inv(A1'*A1);

tempMatrix=eye(num_measurements*2)-A1*inv_A1A1*A1';
D=A2'*tempMatrix*A2;

d=A2'*tempMatrix*b;

lambda=getLambda(D,d,gravity_mag);

I0_g=inv(D-lambda*eye(3))*d;

x1=-inv_A1A1*A1'*A2*I0_g+inv_A1A1*A1'*b;

x=M\b;

%SS=(M*[x1;I0_g]-b)'*(M*[x1;I0_g]-b)

I0_v_I0=x1(num_features*3+1:num_features*3+3,1);

% G_I0_q=getInitQuaternionfromAcc_xfront(I0_g);
% G_I0_R=quatern2rotMat(G_I0_q);

G_I0_R = gram_schmidt(I0_g)';
G_I0_q = rotMat2qRichard(G_I0_R);

for n=1:num_frames

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


validFeatrues = true(num_features, 1); % initial valid Index

for i=1:num_features

    I0_p_f(:,i)=x1(i*3-2:i*3);
    G_p_f(:,i)=G_I0_R*I0_p_f(:,i);

    feat=features(num2str(map_features_ids(i)));

    for cam_id=1:size(feat.timestamps,2)

        for j=1:size(feat.timestamps{cam_id},1)

            camtime=feat.timestamps{cam_id}(j,1);

            n=find_cam_n_from_map_camera_times(camtime,map_camera_times);

            G_Ik_q=x_I_k(1:4,n);

            G_p_Ik=x_I_k(5:7,n);

            G_Ik_R=quatern2rotMat(G_Ik_q);

            G_Ik_Rcam=G_Ik_R*camR;

            G_Ik_Tcam=G_p_Ik+G_Ik_R*camT;

            Vtemp=G_Ik_Rcam'*(G_p_f(:,i)-G_Ik_Tcam);

            if Vtemp(3)<0

                validFeatrues(i,1)=false;
                break;
            end

        end

    end

end
%%

%drawProjection(x_I_k,G_p_f,camR,camT);

G_p_f=G_p_f(:,validFeatrues);

drawProjection(x_I_k,G_p_f,camR,camT);

[features,num_measurements]=features_eliminate_from_validFeatrues(features,validFeatrues,map_features_ids);


%%

[x_I_k_opti,G_p_f_opti]=optimaization_openvins_init_using_features(x_I_k,G_p_f,imuPropagate,features,camK,camD,camR,camT,gravity_mag,map_camera_times,num_measurements);

%%


drawProjection(x_I_k_opti,G_p_f_opti,camR,camT);


%%


frame1=1;
frame2=3;

drawOpticalFlowLK_featrues(imgpyr,features,map_camera_times,cam_id,cam_id,frame1,frame2);

