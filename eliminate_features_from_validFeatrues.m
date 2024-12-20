function [features,num_measurements]=eliminate_features_from_validFeatrues(features,validFeatrues,map_features_ids)

num_measurements=0;

for i = 1:length(map_features_ids)

    id=num2str(map_features_ids(i));

    if validFeatrues(i)==false

        remove(features, id);

    else

        feat=features(id);

        if size(feat.timestamps,2)>0

             for cam_id=1:size(feat.timestamps,2)

                 num_measurements=num_measurements+size(feat.timestamps{cam_id},1);

             end

        end

    end

end


end