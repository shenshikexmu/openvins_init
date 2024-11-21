clear all
clc
addpath('ShanzhaiCV');
addpath('YAMLMatlab');



global matlab_or_octave
matlab_or_octave=0; 

    
file_cam0='../bag/V1_02_medium/mav0/cam0/';

if  (matlab_or_octave ==1)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  matlab   
    datacsv_cam0=readcell([file_cam0,'data.csv']);
    cam0Para = ReadYaml_matlab([file_cam0,'sensor.yaml']);            % cell2mat
    cam0Para.intrinsics=cell2mat(cam0Para.intrinsics);
    cam0Para.distortion_coefficients=cell2mat(cam0Para.distortion_coefficients);
else                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  octave
    pkg load io
    pkg load image
    datacsv_cam0 = csv2cell([file_cam0,'data.csv']);
    cam0Para = readYaml_octave([file_cam0,'sensor.yaml']);
end


global min_px_dist  grid_x grid_y num_features threshold currid



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

max_focallength=max(camK(1,1),camK(2,2));



i=289;



m=10;
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
    ids{n-1}=refine(ids{n-1},mask);
    pts{n}=refine(pts{n},mask);
    ids{n}=ids{n-1};
    
    if n==m
        save(['pts_ids_',num2str(i),'_',num2str(n),'.mat'],'pts','ids');
    end
    
end





for n=1:10
      
    colors = jet(6);
    
    figure;
    imshow([imgpyr{1}{1},imgpyr{n}{1};imgpyr{n}{1},imgpyr{n}{1}]);
    hold on;
    
    s=0;
    for k=1:size(ids{n},1)   
        for j=n:size(ids{1},1)
            if ids{1}(j,1)==ids{n}(k,1)
                
                s=s+1;
                hold on
                plot(pts{1}(j,1),pts{1}(j,2),'r*');
                
                if mod(s,2)==0
                    hold on;
                    plot(pts{n}(k,1),pts{n}(k,2)+size(imgpyr{1}{1},1),'g*');
                    hold on;
                    plot([pts{1}(j,1),pts{n}(k,1)],[pts{1}(j,2),pts{n}(k,2)+size(imgpyr{1}{1},1)],'-', 'Color', colors(mod(s,6)+1,:));
                else
                    hold on;
                    plot(pts{n}(k,1)+size(imgpyr{1}{1},2),pts{n}(k,2),'g*');
                    hold on;
                    plot([pts{1}(j,1),pts{n}(k,1)+size(imgpyr{1}{1},2)],[pts{1}(j,2),pts{n}(k,2)],'-', 'Color', colors(mod(s,6)+1,:));
                end

                hold on;
                quiver(pts{1}(j,1)+size(imgpyr{1}{1},2), pts{1}(j,2)+size(imgpyr{1}{1},1), pts{n}(k,1)-pts{1}(j,1), pts{n}(k,2)-pts{1}(j,2), 'r', 'LineWidth', 2);

            end

        end
        
    end
    
    title(['Optical Flow (Lucas-Kanade) image',num2str(1),' to image',num2str(n)]);
    
    hold off;
    
    a=10;

    
end












