function features=features_p_FinA_from_frame1_frame2(features,map_camera_times,cam_id1,cam_id2,frame1,frame2,R1,T1,R2,T2)


timestamp1=map_camera_times(frame1,1)-map_camera_times(frame1,3);
timestamp2=map_camera_times(frame2,1)-map_camera_times(frame2,3);

ids_all = keys(features);


for i = 1:length(ids_all)

    id = ids_all{i}; % 获取当前键

    feat = features(id);

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

        U1=[feat.uvs_norm{cam_id1}(frame1_timestamps,1);feat.uvs_norm{cam_id1}(frame1_timestamps,2);1];

        U2=[feat.uvs_norm{cam_id2}(frame2_timestamps,1);feat.uvs_norm{cam_id2}(frame2_timestamps,2);1];

        A=[Skew_symmetric(R1*U1);Skew_symmetric(R2*U2)];

        b=[Skew_symmetric(R1*U1)*T1;Skew_symmetric(R2*U2)*T2];

        X=A\b;

        feat.p_FinA=X;

    end

    features(id)=feat;



end



end