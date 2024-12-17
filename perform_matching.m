function mask_out=perform_matching(img0pyr, img1pyr, kpts0, kpts1, id0, id1)
    
    
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

    
    
   
  
 



























 
end
