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

for r=1:1%size(valid_locs,1)

    % Calculate what cell xy value we are in
    grid=valid_locs(r,:);
    x=grid(1)*size_x;
    y=grid(2)*size_y;
    
    % Skip if we are out of bounds
    if (x + size_x > size(img,2) || y + size_y > size(img,1))
        continue
    end
    
    
    
end

 

%
%                      // Calculate where we should be extracting from
%                      cv::Rect img_roi = cv::Rect(x, y, size_x, size_y);
%
%                      // Extract FAST features for this part of the image
%                      std::vector<cv::KeyPoint> pts_new;
%                      cv::FAST(img(img_roi), pts_new, threshold, nonmaxSuppression);
%
%                      // Now lets get the top number from this
%                      std::sort(pts_new.begin(), pts_new.end(), Grider_FAST::compare_response);
%









%                      // Append the "best" ones to our vector
%                      // Note that we need to "correct" the point u,v since we extracted it in a ROI
%                      // So we should append the location of that ROI in the image
%                      for (size_t i = 0; i < (size_t)num_features_grid && i < pts_new.size(); i++) {
%
%                        // Create keypoint
%                        cv::KeyPoint pt_cor = pts_new.at(i);
%                        pt_cor.pt.x += (float)x;
%                        pt_cor.pt.y += (float)y;
%
%                        // Reject if out of bounds (shouldn't be possible...)
%                        if ((int)pt_cor.pt.x < 0 || (int)pt_cor.pt.x > img.cols || (int)pt_cor.pt.y < 0 || (int)pt_cor.pt.y > img.rows)
%                          continue;
%
%                        // Check if it is in the mask region
%                        // NOTE: mask has max value of 255 (white) if it should be removed
%                        if (mask.at<uint8_t>((int)pt_cor.pt.y, (int)pt_cor.pt.x) > 127)
%                          continue;
%                        collection.at(r).push_back(pt_cor);
%                      }
%                    }
%                  }));
%
%    // Combine all the collections into our single vector
%    for (size_t r = 0; r < collection.size(); r++) {
%      pts.insert(pts.end(), collection.at(r).begin(), collection.at(r).end());
%    }
%
%    // Return if no points
%    if (pts.empty())
%      return;
%
%    // Sub-pixel refinement parameters
%    cv::Size win_size = cv::Size(5, 5);
%    cv::Size zero_zone = cv::Size(-1, -1);
%    cv::TermCriteria term_crit = cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 20, 0.001);
%
%    // Get vector of points
%    std::vector<cv::Point2f> pts_refined;
%    for (size_t i = 0; i < pts.size(); i++) {
%      pts_refined.push_back(pts.at(i).pt);
%    }
%
%    // Finally get sub-pixel for all extracted features
%    cv::cornerSubPix(img, pts_refined, win_size, zero_zone, term_crit);
%
%    // Save the refined points!
%    for (size_t i = 0; i < pts.size(); i++) {
%      pts.at(i).pt = pts_refined.at(i);
%    }
%  }
%};
% 
% 
% 
% 
% 
 
 
 
 
 
 
end
