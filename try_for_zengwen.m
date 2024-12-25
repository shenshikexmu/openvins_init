clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');
addpath('repropagate');
addpath('quaternion');

file_cam0='./data/';

global matlab_or_octave
matlab_or_octave=1; 


global min_px_dist  grid_x grid_y num_features threshold currid camK camD win_size pyr_levels

num_pts=800;  % 定位需要提取的特征点数
num_cameras=1;
num_features=round(num_pts/num_cameras);

fast_threshold=20;
threshold=fast_threshold;


camK=[520.9,0,325.1;...
      0,521.0,249.7;...
      0,0,1];

camD=[0,0,0,0];

min_px_dist=10;
grid_x=5; 
grid_y=5;
currid=0;
win_size = [15, 15];
pyr_levels=5;


features=containers.Map();
cam_id=1;

i=1;

m=2;
%table_img_timestamp=zeros(m,2);

for n=1:m
    
    img{n}=rgb2gray(imread([file_cam0,num2str(n),'.png']));
    img{n}=cv_equalizeHist(img{n});
    imgpyr{n}=cv_buildOpticalFlowPyramid(img{n},win_size,pyr_levels);
    timestamp=n;
    %table_img_timestamp(n,:)=[n,timestamp];
    mask{n}=getMask(img{n});
      
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

    end

    for s=1:size(pts{n},1)
        npt_l=undistort_cv(pts{n}(s,1:2)-[1,1], camK,camD);   
        features=update_feature(features, ids{n}(s,1), timestamp, cam_id, pts{n}(s,1), pts{n}(s,2), npt_l(1,1), npt_l(1,2));
    end

end


map_camera_times=[1,1,0;2,2,0];



%%

frame1=1;
frame2=2;

[n,pts1_n,pts2_n,pts1,pts2]=features_intersection_in_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2);

R1=eye(3);
T1=[0;0;0];

[R2_,T2_]=Initial_R_T(pts1_n,pts2_n);

%[R2_,T2_]=optimaization_R_T(R2,T2,pts1_n,pts2_n,pts1,pts2,camK,camD);

features=features_p_FinA_from_frame1_frame2(features,map_camera_times,cam_id,cam_id,frame1,frame2,R1,T1,R2_,T2_);

drawOpticalFlowLK_featrues(imgpyr,features,map_camera_times,cam_id,cam_id,frame1,frame2);

draw_init(features,map_camera_times,R1,T1,R2_,T2_,cam_id,cam_id,frame1,frame2);





























