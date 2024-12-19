function [pts0,ids0]=perform_detection_monocular(img0pyr,mask0,pts0,ids0)
 
 
% @brief Detects new features in the current image
% @param img0pyr              image we will detect features on (first level of pyramid)
% @param mask0                mask which has what ROI we do not want features in
% @param pts0                 vector of currently extracted keypoints in this image
% @param ids0                 vector of feature ids for each currently extracted keypoint
%
% Given an image and its currently extracted features, this will try to add new features if needed.
% Will try to always have the "max_features" being tracked through KLT at each timestep.
% Passed images should already be grayscaled.


global min_px_dist grid_x grid_y num_features threshold currid

if  0 % exist('perform_detection_monocular.mat')==2
    
    load('perform_detection_monocular.mat');
    
else
    
    [rows, cols, ~] = size(img0pyr{1});

    dataType = class(img0pyr{1});  % 确认类型

    size_close = [floor(rows/min_px_dist),floor(cols/min_px_dist)];

    grid_2d_close = cast(zeros(size_close(1), size_close(2), 1),dataType);
    
    size_x=cols/grid_x;
    
    size_y=rows/grid_y;
    
    size_grid = [grid_y,grid_x];

    grid_2d_grid = cast(zeros(size_grid(1),size_grid(2),1),dataType);

    mask0_updated = mask0;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if size(pts0,1)>0

        edge = 10;
        validIdx = true(size(pts0,1), 1); % initial valid Index

        for i = 1:size(pts0,1)
            %Get current left keypoint, check that it is in bounds
            kpt = pts0(i, :);
            x = round(kpt(1));
            y = round(kpt(2));

            if x < edge || x >= cols - edge || y < edge || y >= rows - edge
                validIdx(i) = false;
                continue;
            end
            
            % Calculate mask coordinates for close points
            x_close = floor(kpt(1) / min_px_dist);
            y_close = floor(kpt(2) / min_px_dist);
            
            if x_close < 0 || x_close >= size_close(2) || y_close < 0 || y_close >= size_close(1)
                validIdx(i) = false;
                continue;
            end
            
            % Calculate what grid cell this feature is in
            x_grid = floor(kpt(1) / size_x);
            y_grid = floor(kpt(2) / size_y);
            
            if x_grid < 0 || x_grid >= size_grid(2) || y_grid < 0 || y_grid >= size_grid(1)
                validIdx(i) = false;
                continue;
            end
            
            % Check if this keypoint is near another point
            if grid_2d_close(y_close , x_close ) > 127
                validIdx(i) = false;
                continue;
            end
            
            % Now check if it is in a mask area or not
            % NOTE: mask has max value of 255 (white) if it should be
            if mask0(y , x ) > 127
                validIdx(i) = false;
                continue;
            end
            
            %  Else we are good, move forward to the next point
            grid_2d_close(max(floor(y_close),1) , max(floor(x_close),1) ) = 255;
            if grid_2d_grid(max(floor(y_grid),1) , max(floor(x_grid),1) ) < 255
                grid_2d_grid(max(floor(y_grid),1), max(floor(x_grid),1)) = grid_2d_grid(max(floor(y_grid),1) , max(floor(x_grid),1) ) + 1;
            end
            
            %  Append this to the local mask of the image
            pt1 = [x - min_px_dist, y - min_px_dist];
            pt2 = [x + min_px_dist, y + min_px_dist];
            
            if pt1(1) >= 1 && pt1(2) >= 1 && pt2(1) <= cols && pt2(2) <= rows
                mask0_updated(pt1(2):pt2(2), pt1(1):pt2(1)) = 255;
            end
        end

        pts0 = pts0(validIdx, :);
        ids0 = ids0(validIdx);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    min_feat_percent=0.5;

    num_featsneeded=num_features-size(pts0,1);

    if (num_featsneeded< min(20,min_feat_percent*num_features))
        return
    end


    mask0_grid=cv_resize(mask0, [grid_y,grid_x],0,0, 'INTER_NEAREST');

    num_features_grid= floor(num_features/(grid_x*grid_y))+1;

    num_features_grid_req = max(1, floor(min_feat_percent * num_features_grid));


    valid_locs=[];

    for x=1:size(grid_2d_grid,2)

        for y=1:size(grid_2d_grid,1)
        
            if (grid_2d_grid(y,x)<num_features_grid_req && mask0_grid(y,x)~=255)
             
                valid_locs=[ valid_locs;[x,y]];
             
            end
          
        end

    end

    pts0_ext=perform_griding(img0pyr{1}, mask0_updated, valid_locs, num_features,  grid_x,  grid_y,  threshold, 1) ;

    save('perform_detection_monocular.mat');

end


% Now, reject features that are close a current feature

kpts0_new=[];
pts0_new=[];
for i=1:size(pts0_ext,1)
    
    x_grid=pts0_ext(i,1)/min_px_dist;
    y_grid=pts0_ext(i,2)/min_px_dist;
    
    if (x_grid<0.7 || x_grid> size_close(2) || y_grid <0.7 || y_grid >size_close(1))
        continue;
    end
%    i
%    max(floor(y_grid),1)
%    max(floor(x_grid),1)
    %  See if there is a point at this location
    if (grid_2d_close(max(floor(y_grid),1),max(floor(x_grid),1))>127)
        continue;
    end
    
    
    
    % Else lets add it!
    kpts0_new=[kpts0_new;pts0_ext(i,:)];
    pts0_new=[pts0_new;pts0_ext(i,1:2)];
    grid_2d_close(max(floor(y_grid),1),max(floor(x_grid),1))=255;
   
    
end


% Loop through and record only ones that are valid
% NOTE: if we multi-thread this atomic can cause some randomness due to multiple thread detecting features
% NOTE: this is due to the fact that we select update features based on feat id
% NOTE: thus the order will matter since we try to select oldest (smallest id) to update with
% NOTE: not sure how to remove... maybe a better way?
for i=1:size(pts0_new,1)
    
    % update the uv coordinates
    kpts0_new(i,1:2) = pts0_new(i,1:2);
    % append the new uv coordinate
    pts0=[pts0;kpts0_new(i,:)];
    
    % move id foward and append this new point
    currid=currid+1;
    %temp=currid;
    ids0=[ids0;currid];
   

end

 
  
end
