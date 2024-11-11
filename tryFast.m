%clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  octave
pkg load image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_org=imread('1403715550312143104.png');

[rows, cols] = size(image_org);

image=image_org(rows/2+1:rows,cols/2+1:cols,:);
% 
% image_new=image(1:57,138:246,:);
% 
% imwrite(image_new, 'image_new.png');


threshold=20;

nonmaxSuppression=1;

pts = cv_FAST(image, threshold, nonmaxSuppression);


%pts =[31,9,10];



%    figure
%    imshow(image); hold on;
%    plot(corners(:, 1), corners(:, 2), 'r+'); 


win_size = [5, 5];         % 窗口大小
zero_zone = [-1, -1];      % 禁用区域
term_crit.type = "CV_TERMCRIT_ITER+EPS"; % 终止条件 [最大迭代次数, epsilon]
term_crit.epsilon=0.001;
term_crit.max_iter=20;
   


%  cv::Size win_size = cv::Size(5, 5);
%     cv::Size zero_zone = cv::Size(-1, -1);
%     cv::TermCriteria term_crit = cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 20, 0.001);
% 



pts_refined = cv_cornerSubPix(image, pts(:,1:2), win_size, zero_zone, term_crit);

%%
figure;
imshow(image); hold on;
plot(pts_refined(:,1), pts_refined(:,2), 'g.', 'MarkerSize', 7);
plot(pts(:,1), pts(:,2), 'r.', 'MarkerSize', 7);
%plot(29.57, 11.6, 'b.', 'MarkerSize', 5);
legend('Refined Corners','Initial Corners');
