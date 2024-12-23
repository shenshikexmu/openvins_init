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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  paper <Vins-mono: A robust and versatile monocular visual-inertial state estimator> 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


for n=1:size(map_camera_times,1)

    img{n}=imread([file_cam0,'data/',datacsv_cam0{map_camera_times(n,2),2}]);
    img{n}=cv_equalizeHist(img{n});
    imgpyr{n}=cv_buildOpticalFlowPyramid(img{n},win_size,pyr_levels);

end


for i=1:size(map_camera_times,1)-1

    frame1=i;
    frame2=i+1;

    n=features_intersection_in_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2)

end


%%
frame1=3;
frame2=5;


[n,pts1_n,pts2_n]=features_intersection_in_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2);


R1=eye(3);
T1=[0;0;0];

[R2,T2]=Initial_R_T(pts1_n,pts2_n);

features=features_p_FinA_from_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2,R1,T1,R2,T2);

draw_init(features,R1,T1,R2,T2);

drawOpticalFlowLK_featrues(imgpyr,features,map_camera_times,cam_id,cam_id,frame1,frame2);
































