function features=eliminate_1point_features(features)


allKeys = keys(features);


for i = 1:length(allKeys)

    key = allKeys{i}; % 获取当前键

    feat = features(key);

    sum_point=0;

    if size(feat.timestamps,2)>0

        for cam_id=1:size(feat.timestamps,2)

            sum_point=sum_point+size(feat.timestamps{cam_id},1);

        end

    end

    if sum_point<2
        remove(features, key);
        
    end




end


end