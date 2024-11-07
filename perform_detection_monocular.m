function [pts0,ids0]=perform_detection_monocular(img0pyr,mask0,pts0,ids0)
  
  
global min_px_dist grid_x grid_y num_features threshold

[rows, cols, ~] = size(img0pyr{1});

dataType = class(img0pyr{1});  % 确认类型

%size_close=[int(rows/min_px_dist),int(cols/min_px_dist)];

grid_2d_close = cast(zeros(floor(rows/min_px_dist), floor(cols/min_px_dist), 1),dataType);

%size_grid=[grid_y,grid_x];

grid_2d_grid= cast(zeros(grid_y,grid_x,1),dataType);

mask0_updated=mask0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 没有写
if size(pts0,1)>0

    for i=1:size(pts0,1)
    
    
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_feat_percent=0.5;

num_featsneeded=num_features-size(pts0,1);

if (num_featsneeded< min(20,min_feat_percent*num_features))
    return
end


mask0_grid=cv_resize(mask0, [grid_y,grid_x],0,0, 'INTER_NEAREST');

num_features_grid= floor(num_features/(grid_x*grid_y))+1;

num_features_grid_req = max(1, floor(min_feat_percent * num_features_grid))




valid_locs=[];

for x=1:size(grid_2d_grid,2)

    for y=1:size(grid_2d_grid,1)
    
        if (grid_2d_grid(y,x)<num_features_grid_req && mask0_grid(y,x)~=255)
         
            valid_locs=[ valid_locs;[x,y]];
         
        end
    
    
    end

end




pts0_ext=perform_griding(img0pyr{1}, mask0_updated, valid_locs, num_features,  grid_x,  grid_y,  threshold, 1) ;

















  
  
  
  
  
  
  
end
