clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');


global matlab_or_octave
matlab_or_octave=0; 

    
file_cam0='../bag/V1_02_medium/mav0/cam0/';


if  (matlab_or_octave ==1)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  matlab   
    cam0Para = ReadYaml_matlab([file_cam0,'sensor.yaml']);            % cell2mat
    cam0Para.intrinsics=cell2mat(cam0Para.intrinsics);
    cam0Para.distortion_coefficients=cell2mat(cam0Para.distortion_coefficients);
    cam0Para.T_BS.data=reshape(cell2mat(cam0Para.T_BS.data),4,4)';
else                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  octave
    pkg load io
    pkg load image
   
    cam0Para = readYaml_octave([file_cam0,'sensor.yaml']);
    cam0Para.T_BS.data=reshape(cam0Para.data,4,4)';
    
end




camK=[cam0Para.intrinsics(1),0,cam0Para.intrinsics(3);...
      0,cam0Para.intrinsics(2),cam0Para.intrinsics(4);...
      0,0,1];

camD=[cam0Para.distortion_coefficients(1),cam0Para.distortion_coefficients(2),cam0Para.distortion_coefficients(3),cam0Para.distortion_coefficients(4)];   


uv_norm=[0.1,0.2];
delta=0.00001;

[uv_dist,H_dz_dzn] = distort_cv(uv_norm, camK,camD);


for i=1:2
    
    
    uv_norm_temp=uv_norm;
    
    uv_norm_temp(i)=uv_norm_temp(i)+delta;
    
    [uv_dist_temp,H_dz_dzn] = distort_cv(uv_norm_temp, camK,camD);
    
    J(:,i)=(uv_dist_temp'-uv_dist')/delta;
    
    
    
end



