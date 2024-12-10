clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');

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

camK=[cam0Para.intrinsics(1),0,cam0Para.intrinsics(3);...
      0,cam0Para.intrinsics(2),cam0Para.intrinsics(4);...
      0,0,1];

camD=[cam0Para.distortion_coefficients(1),cam0Para.distortion_coefficients(2),cam0Para.distortion_coefficients(3),cam0Para.distortion_coefficients(4)];

camR=cam0Para.T_BS.data(1:3,1:3);

camT=cam0Para.T_BS.data(1:3,4);



features=get_features_from_txt('./data/feature.txt');

map_camera_times = [];
count_valid_features=0;

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
            
            count_valid_features=count_valid_features+size(feat.uvs_norm{cam_id},1);
            
        end
    end
    
    
end

map_camera_times = unique(map_camera_times);


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
