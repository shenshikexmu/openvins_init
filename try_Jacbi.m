%clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');
addpath('quaternion');
addpath('repropagate');

global matlab_or_octave
matlab_or_octave=1; 

    
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


camR=cam0Para.T_BS.data(1:3,1:3);

camT=cam0Para.T_BS.data(1:3,4);




uv_norm=[0.1,0.2];
delta=0.00001;

[uv_dist,H_dz_dzn] = distort_cv(uv_norm, camK,camD);


for i=1:2
    
    
    uv_norm_temp=uv_norm;
    
    uv_norm_temp(i)=uv_norm_temp(i)+delta;
    
    [uv_dist_temp,H_dz_dzn] = distort_cv(uv_norm_temp, camK,camD);
    
    J(:,i)=(uv_dist_temp'-uv_dist')/delta;
    
    
    
end


theta=randn(3,1);
theta0=randn(3,1);%[0;0;0];%randn(3,1);
a=randn(3,1);
%a=a/norm(a);

R0=angleAxisToRotationMatrix(theta0);

for i=1:3
    
    theta_temp=theta;
    theta_temp(i)=theta_temp(i)+delta;
    
    JRa(:,i)=(R0*angleAxisToRotationMatrix(theta_temp)*a-R0*angleAxisToRotationMatrix(theta)*a)/delta;
    
    JRtransa(:,i)=(R0'*angleAxisToRotationMatrix(theta_temp)'*a-R0'*angleAxisToRotationMatrix(theta)'*a)/delta;
    
end

%JRa
%
%dRa_dtheta=R0*D_Ra_D_theta(theta,a)
%
%JRtransa
%
%dRtans_a_dtheta=R0'*D_Rtrans_a_D_theta(theta,a)


theta=randn(3,1);
theta0=randn(3,1);%[0;0;0];%randn(3,1);
q=angleAxis2Quaternion(theta);
q0=angleAxis2Quaternion(theta0);

a=randn(3,1);
%a=a/norm(a);

R0=quatern2rotMat(q0);

deltaX=[0;0;0];

for i=1:3
    
    deltaX_temp=deltaX;

    deltaX_temp(i)=deltaX_temp(i)+delta;



    q_temp=quaternProd(q, utility_deltaQ_VINS_Mono(deltaX_temp) );
    
    
    JRa(:,i)=(R0*quatern2rotMat(q_temp)*a-R0*quatern2rotMat(q)*a)/delta;
    
    JRtransa(:,i)=(R0*quatern2rotMat(q_temp)'*a-R0*quatern2rotMat(q)'*a)/delta;
    
end

% JRa
% 
% dRa_ddelta_theta=R0*D_Ra_D_delta_theta(q,a)
% 
% JRtransa
% 
% dRtans_a_ddelta_theta=R0*D_Rtransa_D_delta_theta(q,a)







G_I_p=[0;0;0];

G_I_angleAxis=[0.051814;-1.9161;-0.073786];

G_I_q=angleAxis2Quaternion(G_I_angleAxis);

G_p_f_k=[-0.042518;-0.084215;-0.050919];

pts_temp=[144.09,	316.37];

[E,H_dz_dG_I_p,H_dz_dG_I_angleAxis,H_dz_dG_p_f_k]=evaluate_Reprojection_try(G_I_p,G_I_q,G_p_f_k,camK,camD,camR,camT,pts_temp);

deltaX=[0;0;0];

for i=1:3

    deltaX_temp=deltaX;

    deltaX_temp(i)=deltaX_temp(i)+delta;



    G_I_q_temp=quaternProd(G_I_q, utility_deltaQ_VINS_Mono(deltaX_temp) );
    
    [E_temp,H_dz_dG_I_p,H_dz_dG_I_angleAxis,H_dz_dG_p_f_k]=evaluate_Reprojection_try(G_I_p,G_I_q_temp,G_p_f_k,camK,camD,camR,camT,pts_temp);

    JG_I_angleAxis(:,i)=(E_temp-E)/delta;
    
end


JG_I_angleAxis
 

H_dz_dG_I_angleAxis








