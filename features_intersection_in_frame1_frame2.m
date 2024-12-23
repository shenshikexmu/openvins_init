function [n,pts1_n,pts2_n]=features_intersection_in_frame1_frame2(features,map_camera_times,cam_id1,cam_id2,frame1,frame2)


timestamp1=map_camera_times(frame1,1)-map_camera_times(frame1,3);
timestamp2=map_camera_times(frame2,1)-map_camera_times(frame2,3);

allKeys = keys(features);
n=0;

pts1_n=[];
pts2_n=[];


for i = 1:length(allKeys)

    key = allKeys{i}; % 获取当前键

    feat = features(key);

    frame1_timestamps=0;

    if size(feat.timestamps{cam_id1},1)>1

        for j=1:size(feat.timestamps{cam_id1},1)

            if timestamp1==feat.timestamps{cam_id1}(j)
                
                frame1_timestamps=j;

            end

        end

    end

    frame2_timestamps=0;

    if size(feat.timestamps{cam_id2},1)>1

        for j=1:size(feat.timestamps{cam_id2},1)

            if timestamp2==feat.timestamps{cam_id2}(j)
                
                frame2_timestamps=j;

            end

        end

    end

    if frame1_timestamps~=0 && frame2_timestamps~=0

        n=n+1;

        pts1_n=[pts1_n;feat.uvs_norm{cam_id1}(frame1_timestamps,:)];

        pts2_n=[pts2_n;feat.uvs_norm{cam_id2}(frame2_timestamps,:)];

    end

end



end