function n=find_cam_n_from_map_camera_times(camtime,map_camera_times)
    
min_gap=Inf;

for i=1:size(map_camera_times,1)
    
    gap=abs(map_camera_times(i,1)-camtime);
    
    if gap<min_gap
        
        min_gap=gap;
        
        n=i;
    end

end
 
end
