function [pts1, status, err] = cv_calcOpticalFlowPyrLK(imgpyr0, imgpyr1, pts0, pts1,win_size, pyr_levels, criteria, cv_OPTFLOW_USE_INITIAL_FLOW)

max_iters=criteria.max_iters;
eps=criteria.epsilon;

if (cv_OPTFLOW_USE_INITIAL_FLOW==1)
    pts1_=pts1;
else
    pts1_=pts0;
end

for n=pyr_levels:-1:1
    
    I0_now=imgpyr0{n};
    I1_now=imgpyr1{n};
    
    scale=2^(n-1);
    
    pts0_now=pts0/scale;
    pts1_now=pts1/scale;
    

%    figure
%    imshow([I0_now;I1_now]);
%    hold on;
%    plot(pts0_now(:,1),pts0_now(:,2),'r.');
%     
    
    [pts1] = iterOpticalflow(I0_now,I1_now,pts0_now,pts1_now,win_size,eps,max_iters);
    
    pts1=pts1*scale;
    
    a=10;
    
end


for i=1:size(pts0,1)
    
    figure
    imshow([imgpyr0{1},imgpyr1{1}]);
    hold on;
    plot(pts0(i,1),pts0(i,2),'r.', 'MarkerSize', 7,pts1(i,1)+size(imgpyr0{1},2),pts1(i,2),'r.', 'MarkerSize', 7);
    
    
    a=10;
    
    
end

a=10;
                               
                               
end