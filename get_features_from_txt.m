function features=get_features_from_txt(file_name)
    
    
feature_txt=importdata(file_name);

features=containers.Map();

for i=1:size(feature_txt,1)
    
    featid=feature_txt(i,1)+1;
    cam_id=feature_txt(i,2)+1;
    timestamps=feature_txt(i,3);
    uvs_norm=feature_txt(i,4:5);
    uvs=feature_txt(i,6:7);
    
    feat=[];
    if ~isKey(features,num2str(featid))

        feat.timestamps=[];    
        feat.uvs_norm{cam_id}=[];     
        feat.uvs{cam_id}=[];    
        feat.featid=featid;
        
        features(num2str(featid))=feat;
        
    end
    
    feat=features(num2str(featid));
        
    feat.timestamps=[feat.timestamps;timestamps];
    
    feat.uvs_norm{cam_id}=[feat.uvs_norm{cam_id};uvs_norm];
    
    feat.uvs{cam_id}=[feat.uvs{cam_id};uvs];
        
    features(num2str(featid))=feat;

end
end
