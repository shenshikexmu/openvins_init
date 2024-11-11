%clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  matlab
file_cam0='../bag/V1_02_medium/mav0/cam0/';
datacsv_cam0=readcell([file_cam0,'data.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% octave
% pkg load io
% pkg load image
% file_cam0 ='../bag/V1_02_medium/mav0/cam0/';
% datacsv_cam0 = csv2cell([file_cam0,'data.csv']);



global min_px_dist grid_x grid_y num_features threshold currid


min_px_dist=10;
grid_x=5;
grid_y=5;

pyr_levels=5;


num_pts=200;  % 定位需要提取的特征点数
num_cameras=1;

num_features=round(num_pts/num_cameras);

fast_threshold=20;

threshold=fast_threshold;

currid=0;

win_size = [15, 15];


i=700;

img0=imread([file_cam0,'data/',datacsv_cam0{i,2}]);
    
mask=getMask(img0);

imgpyr0=buildOpticalFlowPyramid(img0,win_size,pyr_levels);


pts0=[];
ids0=[];

if (exist('pts0_ids0.mat')==2)
    load('pts0_ids0.mat');
else
    [pts0,ids0]=perform_detection_monocular(imgpyr0,mask,pts0,ids0);
    save('pts0_ids0.mat','pts0','ids0');
end

% figure 
% imshow(imgpyr0{1});hold on;
% plot(pts0(:,1), pts0(:,2), 'r.', 'MarkerSize', 7); 


img1=imread([file_cam0,'data/',datacsv_cam0{i+0,2}]);
imgpyr1=buildOpticalFlowPyramid(img1,win_size,pyr_levels);


    %cv::TermCriteria term_crit = cv::TermCriteria(cv::TermCriteria::COUNT | cv::TermCriteria::EPS, 30, 0.01);
% criteria=[30,0.01];
% 
%  %cv::OPTFLOW_USE_INITIAL_FLOW
%  derivLambda=[];
%  flags=0;
% 
% %[nextPts, status, err] = calcOpticalFlowPyrLK_doubao(imgpyr0, imgpyr1, pts0, win_size, pyr_levels, criteria, derivLambda, flags)


% 调用Lucas-Kanade光流算法
window_size = 5; % 设置窗口大小
[u, v] = lucas_kanade_optical_flow_points(imgpyr0{1}, imgpyr1{1}, pts0, window_size);

colors = jet(6);
figure;
imshow([imgpyr0{1};imgpyr1{1}]);
hold on;
for i=1:size(pts0,1)
    hold on
    plot([pts0(i,1),pts0(i,1)+u(i)],[pts0(i,2),pts0(i,2)+480+v(i)],'-', 'Color', colors(mod(i,6)+1,:));

end
%quiver(pts0(:,1), pts0(:,2), u+pts0(:,1)+752, v, 'r');
title('Optical Flow (Lucas-Kanade at Feature Points)');



%for i=800:800
%
%    img=imread([file_cam0,'data/',datacsv_cam0{i,2}]);
%    
%    mask=getMask(img);
%    
%    imgpyr=buildOpticalFlowPyramid(img,[15,15],1);
%    
%    pts0=[];
%    ids0=[];
%    
%    if (exist('pts0_ids0.mat')==2)
%        load('pts0_ids0.mat');
%    else
%        [pts0,ids0]=perform_detection_monocular(imgpyr,mask,pts0,ids0);
%        save('pts0_ids0.mat','pts0','ids0');
%    end
%
%   
%    figure 
%
%    imshow(imgpyr{1});
%
%    hold on;
%    plot(pts0(:,1), pts0(:,2), 'r.', 'MarkerSize', 7); 
%
%
%end

%%


%% 读取图像帧
%I1 = imread([file_cam0,'data/',datacsv_cam0{1000,2}]);
%I2 = imread([file_cam0,'data/',datacsv_cam0{1005,2}]);
%
%% % 使用 Lucas-Kanade 算法计算光流
%% opticFlow = opticalFlowLK('NoiseThreshold',0.01);
%% flow1 = estimateFlow(opticFlow, I1);
%% flow2 = estimateFlow(opticFlow, I2);
%% 
%% % % 使用 Horn-Schunck 算法计算光流
%% % opticFlowHS = opticalFlowHS('Smoothness',1);
%% % flowHS1 = estimateFlow(opticFlowHS, I1);
%% % flowHS2 = estimateFlow(opticFlowHS, I2);
%% 
%% % 显示光流
%% figure;
%% imshow(I1);
%% hold on;
%% plot(flow1, 'DecimationFactor', [5 5], 'ScaleFactor', 10);
%% 
%% 
%
%
%
%% 使用示例
%% I1 = rgb2gray(imread('frame1.png')); % 读取第一帧图像
%% I2 = rgb2gray(imread('frame2.png')); % 读取第二帧图像
%window_size = 5; % 设置卷积窗口的大小
%
%% 调用光流计算函数
%[u, v] = lucas_kanade_optical_flow(I1, I2, window_size);
%
%% 显示光流
%figure;
%imshow(I1); % 显示第一帧图像
%hold on;
%quiver(u, v, 5, 'r'); % 绘制光流，箭头大小可调整
%title('Lucas-Kanade Optical Flow');

