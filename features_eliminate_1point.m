function features=features_eliminate_1point(features)


all_ids = keys(features);


for i = 1:length(all_ids)

    id = all_ids{i}; % 获取当前键

    feat = features(id);

    sum_point=0;

    if size(feat.timestamps,2)>0

        for cam_id=1:size(feat.timestamps,2)

            sum_point=sum_point+size(feat.timestamps{cam_id},1);

        end

    end

    if sum_point<2
        
        remove(features, id);
        
    end

end


end