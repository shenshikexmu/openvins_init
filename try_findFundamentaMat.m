clear all
clc
addpath('ShanzhaiCV');




global matlab_or_octave=0 
global min_px_dist=10  grid_x=5 grid_y=5 num_features threshold currid=0 

    

if  (matlab_or_octave ==1)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  matlab
    file_cam0='../bag/V1_02_medium/mav0/cam0/';
    datacsv_cam0=readcell([file_cam0,'data.csv']);
else                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  octave
    pkg load io
    pkg load image
    file_cam0 ='../bag/V1_02_medium/mav0/cam0/';
    datacsv_cam0 = csv2cell([file_cam0,'data.csv']);
end

cam0Para = readYaml([file_cam0,'sensor.yaml']);



pyr_levels=5;

num_pts=200;  % 定位需要提取的特征点数
num_cameras=1;

num_features=round(num_pts/num_cameras);

fast_threshold=20;

threshold=fast_threshold;

win_size = [15, 15];


i=539;

img0=imread([file_cam0,'data/',datacsv_cam0{i,2}]);
img0=cv_equalizeHist(img0);
imgpyr0=cv_buildOpticalFlowPyramid(img0,win_size,pyr_levels);


img1=imread([file_cam0,'data/',datacsv_cam0{i+1,2}]);
img1=cv_equalizeHist(img1);
imgpyr1=cv_buildOpticalFlowPyramid(img1,win_size,pyr_levels);


mask=getMask(img0);


pts0=[];
ids0=[];

if (exist(['pts0_ids0_pts1_',num2str(i),'.mat'])==2)
    load(['pts0_ids0_pts1_',num2str(i),'.mat']);
else
    [pts0,ids0]=perform_detection_monocular(imgpyr0,mask,pts0,ids0);
    criteria.max_iters=30;
    criteria.epsilon=0.01;
    pts1=pts0;
    cv_OPTFLOW_USE_INITIAL_FLOW=1;
    [pts1, status, err] = cv_calcOpticalFlowPyrLK(imgpyr0, imgpyr1, pts0, pts1, win_size, pyr_levels, criteria, cv_OPTFLOW_USE_INITIAL_FLOW);
    save(['pts0_ids0_pts1_',num2str(i),'.mat'],'pts0','ids0','pts1','status','err');
end



%mask=cv_findFundamentalMat(points1, points2, method, ransacReprojThreshold ,confidence )























%%
figure;
imshow([imgpyr0{1},imgpyr1{1};imgpyr1{1},imgpyr1{1}]);
hold on;
plot(pts0(:,1),pts0(:,2),'r*');



colors = jet(6);

n=0;
for i=1:size(pts0,1)
    
    if status(i,1)==1
        n=n+1;
        hold on;
        plot(pts0(i,1),pts0(i,2),'r*');
        if mod(n,2)==0
            hold on;
            plot(pts1(i,1),pts1(i,2)+size(imgpyr0{1},1),'g*');
            hold on;
            plot([pts0(i,1),pts1(i,1)],[pts0(i,2),pts1(i,2)+size(imgpyr0{1},1)],'-', 'Color', colors(mod(i,6)+1,:));
        else
            hold on;
            plot(pts1(i,1)+size(imgpyr0{1},2),pts1(i,2),'g*');
            hold on;
            plot([pts0(i,1),pts1(i,1)+size(imgpyr0{1},2)],[pts0(i,2),pts1(i,2)],'-', 'Color', colors(mod(i,6)+1,:));
        end

        hold on;
        quiver(pts0(i,1)+size(imgpyr0{1},2), pts0(i,2)+size(imgpyr0{1},1), pts1(i,1)-pts0(i,1), pts1(i,2)-pts0(i,2), 'r', 'LineWidth', 2);
        
        
    end

end

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



