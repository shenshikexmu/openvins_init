function pts=perform_griding(img, mask, valid_locs, num_features,  grid_x,  grid_y,  threshold, nonmaxSuppression) 
 
 
% @brief This function will perform grid extraction using FAST.
% @param img Image                 we will do FAST extraction on
% @param mask                      Region of the image we do not want to extract features in (255 = do not detect features)
% @param valid_locs                Valid 2d grid locations we will extract in (instead of the whole image)
% @param pts                       vector of extracted points we will return
% @param num_features              max number of features we want to extract
% @param grid_x                    size of grid in the x-direction / u-direction
% @param grid_y                    size of grid in the y-direction / v-direction
% @param threshold                 FAST threshold paramter (10 is a good value normally)
% @param nonmaxSuppression         if FAST should perform non-max suppression (true normally)
%
% Given a specified grid size, this will try to extract fast features from each grid.
% It will then return the best from each grid in the return vector.
% 
 
% Return if there is nothing to extract
if (size(valid_locs,1)==0)
    return
end
% We want to have equally distributed features
% NOTE: If we have more grids than number of total points, we calc the biggest grid we can do
% NOTE: Thus if we extract 1 point per grid we have
% NOTE:    -> 1 = num_features / (grid_x * grid_y)
% NOTE:    -> grid_x = ratio * grid_y (keep the original grid ratio)
% NOTE:    -> grid_y = sqrt(num_features / ratio)

if (num_features < grid_x * grid_y) 
    ratio = grid_x / grid_y;
    grid_y = ceil(sqrt(num_features / ratio));
    grid_x = ceil(grid_y * ratio);
end

num_features_grid = num_features / (grid_x * grid_y) + 1;

% Calculate the size our extraction boxes should be
size_x = size(img,2) / grid_x;
size_y = size(img,1) / grid_y;


% Parallelize our 2d grid extraction!!
%collection(valid_locs.size());

collection=[];

for r=1:size(valid_locs,1)

    % Calculate what cell xy value we are in
    grid=valid_locs(r,:);
    x=floor((grid(1)-1)*size_x+1);
    y=floor((grid(2)-1)*size_y+1);
    
    % Skip if we are out of bounds
    if (x + size_x -1 > size(img,2) || y + size_y -1 > size(img,1))
%        x + size_x -1
%        y + size_y -1
        continue;
    end
    
    img_roi=img(y:y+floor(size_y)-1,x:x+floor(size_x)-1);
    
    % cv_FAST(img(y:y+floor(size_y),x:x+floor(size_x)), pts_new, threshold, nonmaxSuppression);
    pts_new = cv_FAST(img_roi, threshold, nonmaxSuppression);

%    figure
%    imshow(img_roi); hold on;
%    plot(pts_new(:, 1), pts_new(:, 2), 'r+'); 

    % Now lets get the top number from this
    pts_new = sortrows(pts_new, -3); 
    
    % Append the "best" ones to our vector
    % Note that we need to "correct" the point u,v since we extracted it in a ROI
    % So we should append the location of that ROI in the image
    collection_tmp=[];
    for i=1:min(num_features_grid,size(pts_new,1))
        %Create keypoint
        pt_cor_pt_x=pts_new(i,1)+x-1;
        pt_cor_pt_y=pts_new(i,2)+y-1;

        % Reject if out of bounds (shouldn't be possible...)

        if ( pt_cor_pt_x < 0 || pt_cor_pt_x > size(img,2) || pt_cor_pt_y < 0 || pt_cor_pt_y> size(img,1))
            continue;
        end

        % Check if it is in the mask region
        % NOTE: mask has max value of 255 (white) if it should be removed
        if (mask(pt_cor_pt_y,pt_cor_pt_x) >127)
            continue;
        end

        collection_tmp=[collection_tmp;pt_cor_pt_x,pt_cor_pt_y,pts_new(i,3)];

    end

    collection{r}=collection_tmp;
    
end

% Combine all the collections into our single vector
pts=[];
for r=1:length(collection)
    
    pts=[pts;collection{r}];

end

% figure
% imshow(img); hold on;
% plot(pts(:, 1), pts(:, 2), 'r+'); 

%  Return if no points
if (size(pts,1)==0)
    return;
end
 
win_size = [5, 5];         % 窗口大小
zero_zone = [-1, -1];      % 禁用区域
term_crit.type = "CV_TERMCRIT_ITER+EPS"; % 终止条件 [最大迭代次数, epsilon]
term_crit.epsilon=0.001;
term_crit.max_iter=20;

% Finally get sub-pixel for all extracted features
pts_refined = cv_cornerSubPix(img, pts(:,1:2), win_size, zero_zone, term_crit);


% figure;
% imshow(img); hold on;
% plot(pts(:,1), pts(:,2), 'r.', 'MarkerSize', 7);
% plot(pts_refined(:,1), pts_refined(:,2), 'g.', 'MarkerSize', 7);
% legend('Initial Corners', 'Refined Corners');

% Save the refined points!
pts(:,1:2)=pts_refined;

 
 
end
