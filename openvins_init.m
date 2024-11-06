clear all
clc


file_cam0='../bag/V1_02_medium/mav0/cam0/';

datacsv_cam0=readcell([file_cam0,'data.csv']);




for i=1000:1000

    image=imread([file_cam0,'data/',datacsv_cam0{i,2}]);

    figure
    imshow(image);

    threshold = 20;  % 设置角点检测的阈值
    keypoints = fast_corner_detection(image, threshold);  % 调用FAST角点检测

    

    hold on;
    plot(keypoints(:,1), keypoints(:,2), 'ro');  % 在图像上绘制角点

    figure

    corners=detectFASTFeatures(image);

    imshow(image); hold on;
    plot(corners);


end

%%


% 读取图像帧
I1 = imread([file_cam0,'data/',datacsv_cam0{1000,2}]);
I2 = imread([file_cam0,'data/',datacsv_cam0{1005,2}]);

% % 使用 Lucas-Kanade 算法计算光流
% opticFlow = opticalFlowLK('NoiseThreshold',0.01);
% flow1 = estimateFlow(opticFlow, I1);
% flow2 = estimateFlow(opticFlow, I2);
% 
% % % 使用 Horn-Schunck 算法计算光流
% % opticFlowHS = opticalFlowHS('Smoothness',1);
% % flowHS1 = estimateFlow(opticFlowHS, I1);
% % flowHS2 = estimateFlow(opticFlowHS, I2);
% 
% % 显示光流
% figure;
% imshow(I1);
% hold on;
% plot(flow1, 'DecimationFactor', [5 5], 'ScaleFactor', 10);
% 
% 



% 使用示例
% I1 = rgb2gray(imread('frame1.png')); % 读取第一帧图像
% I2 = rgb2gray(imread('frame2.png')); % 读取第二帧图像
window_size = 5; % 设置卷积窗口的大小

% 调用光流计算函数
[u, v] = lucas_kanade_optical_flow(I1, I2, window_size);

% 显示光流
figure;
imshow(I1); % 显示第一帧图像
hold on;
quiver(u, v, 5, 'r'); % 绘制光流，箭头大小可调整
title('Lucas-Kanade Optical Flow');

