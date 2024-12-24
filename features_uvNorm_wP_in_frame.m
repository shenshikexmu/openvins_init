function [n,pts_n,pts,wP]=features_uvNorm_wP_in_frame(features,map_camera_times,cam_id,frame)


timestamp=map_camera_times(frame,1)-map_camera_times(frame,3);

allKeys = keys(features);
n=0;

pts_n=[];
pts=[];
wP=[];

for i = 1:length(allKeys)

    key = allKeys{i}; % 获取当前键

    feat = features(key);

    if size(feat.timestamps{cam_id},1)>1

        if ~isempty(feat.p_FinA)

            for j=1:size(feat.timestamps{cam_id},1)
    
                if timestamp==feat.timestamps{cam_id}(j)
    
                    n=n+1;
    
                    pts_n=[pts_n;feat.uvs_norm{cam_id}(j,:)];

                    pts=[pts;feat.uvs{cam_id}(j,:)];

                    wP=[wP;feat.p_FinA'];
    
                end
    
            end
        end

    end

end



end