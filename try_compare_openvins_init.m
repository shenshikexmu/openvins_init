% clear all
% clc
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

camK=[cam0Para.intrinsics(1),0,cam0Para.intrinsics(3);...
      0,cam0Para.intrinsics(2),cam0Para.intrinsics(4);...
      0,0,1];

camD=[cam0Para.distortion_coefficients(1),cam0Para.distortion_coefficients(2),cam0Para.distortion_coefficients(3),cam0Para.distortion_coefficients(4)];

camR=cam0Para.T_BS.data(1:3,1:3);

camT=cam0Para.T_BS.data(1:3,4);

gravity_mag= 9.81;

global ACC_N GYR_N ACC_W GYR_W 


ACC_N=2.0000e-3;                     
GYR_N=1.6968e-04;
ACC_W=3.0000e-3;
GYR_W=1.9393e-05;

features=get_features_from_txt('./data/features.txt');

map_camera_times = [];
count_valid_features=0;
num_measurements=0;

allKeys = keys(features);

for i = 1:length(allKeys)
    key = allKeys{i}; % 获取当前键

    feat = features(key);
    
    map_camera_ids(i,1)=feat.featid;
    
    if ~isempty(feat.timestamps)
        if ~ismember(feat.timestamps,map_camera_times)
            map_camera_times = [map_camera_times; feat.timestamps];
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
    
    for cam_id = 1:length(feat.uvs_norm)
        uvs_norm=feat.uvs_norm{cam_id};
        uvs=feat.uvs{cam_id};
        
        timestamps=feat.timestamps;
        
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


for i=1:num_features

    I0_p_f(:,i)=x1(i*3-2:i*3);
    G_p_f(:,i)=G_I0_R*I0_p_f(:,i);

end

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

drawProjection(x_I_k,G_p_f,camR,camT);



[x_I_k_opti,G_p_f_opti]=optimaization_openvins_init_using_features(x_I_k,G_p_f,imuPropagate,features,camK,camD,camR,camT,gravity_mag,map_camera_times,num_measurements);


drawProjection(x_I_k_opti,G_p_f_opti,camR,camT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isCompare=0;

if isCompare

    M=importdata('./data/A.txt');

    b=importdata('./data/b.txt');

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



    I0_v_I0=x1(num_features*3+1:num_features*3+3,1);

    G_I0_q=getInitQuaternionfromAcc_xfront(I0_g);
    G_I0_R=quatern2rotMat(G_I0_q);

    for i=1:num_features

        I0_p_f(:,i)=x1(i*3-2:i*3);

        G_p_f(:,i)=G_I0_R*I0_p_f(:,i);

    end



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

    drawProjection(x_I_k,G_p_f,camR,camT);

end



%
%
%feature(num2str(5))=10;
%
%feature('8')=20;
%
%a=10
%
%for i=1:size(feature_txt,1)
%    
%    feature(feature_txt(i,1))=
%    
%    
%    if 
%    
%    
%    
%   %    
%    
%    
