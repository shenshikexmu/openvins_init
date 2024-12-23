function features=update_feature(features, id, timestamp, cam_id, u, v, u_n, v_n)
    
feat=[];
if ~isKey(features,num2str(id))

    feat.timestamps{cam_id}=[];    
    feat.uvs_norm{cam_id}=[];     
    feat.uvs{cam_id}=[]; 
    feat.p_FinA=[];
    feat.p_FinG=[];
    feat.featid=id;
    
    features(num2str(id))=feat;
    
end

feat=features(num2str(id));

if cam_id>size(feat.timestamps,2)

    feat.timestamps{cam_id}=[];    
    feat.uvs_norm{cam_id}=[];     
    feat.uvs{cam_id}=[]; 

end

    
feat.timestamps{cam_id}=[feat.timestamps{cam_id};timestamp];

feat.uvs_norm{cam_id}=[feat.uvs_norm{cam_id};[u_n,v_n]];

feat.uvs{cam_id}=[feat.uvs{cam_id};[u,v]];
    
features(num2str(id))=feat;
    
end
