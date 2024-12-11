function [n_frame,timestamps_new, min_gap]=find_image_frame_corresponding_timestamps(datacsv_cam0,timestamps)
    
min_gap=Inf;
    
for i=2:size(datacsv_cam0,1)
    
    gap=abs(datacsv_cam0{i,1}*10e-10-timestamps);
    
    
    if gap<min_gap
        
        min_gap=gap;
        
        n_frame=i;
        
    end
    
    
    
end
    
timestamps_new =  datacsv_cam0{n_frame,1}*10e-10  ;
  
    
end
