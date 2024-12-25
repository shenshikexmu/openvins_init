function  [kpts0, kpts1, mask_out]=perform_matching(img0pyr, img1pyr, kpts0, kpts1, id0, id1)
    
    
    global camK camD win_size pyr_levels
    
    
    if size( kpts0,1)~=size( kpts1,1)
        mask_out=[];
        disp('size( kpts0,1)~=size( kpts1,1)');
        return;
        
    end
    
    if size( kpts0,1)==0
        mask_out=[];
        return;
    end
    
    pts0=kpts0;
    pts1=kpts1;
    
    if size( kpts0,1)<10
        mask_out=zeros(size( kpts0,1),1);
        return;
    end

    term_crit.max_iters=30;
    term_crit.epsilon=0.01;
    cv_OPTFLOW_USE_INITIAL_FLOW=1;
    
 
    [pts1, mask_klt, error] = cv_calcOpticalFlowPyrLK(img0pyr, img1pyr, pts0, pts1, win_size, pyr_levels, term_crit, cv_OPTFLOW_USE_INITIAL_FLOW);

    pts0_n=zeros(size(pts0,1),2);
    pts1_n=zeros(size(pts1,1),2);

    for i=1:size(pts0,1)
        pts0_n(i,:)=undistort_cv(pts0(i,1:2)-[1,1], camK,camD);
        pts1_n(i,:)=undistort_cv(pts1(i,1:2)-[1,1], camK,camD);
    end
    
    max_focallength=max(camK(1,1),camK(2,2));
    
    [mask_rsc]=cv_findFundamentalMat(pts0_n, pts1_n, 'cv_FM_RANSAC', 1/max_focallength ,0.9999);
    
    
    mask_out=floor((mask_klt+mask_rsc)/2);
    
    kpts0=pts0;
    
    kpts1=pts1;
    
end
