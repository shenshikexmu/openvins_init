function features=update_feature(features, id, timestamp, cam_id, u, v, u_n, v_n)
    
    feat=[];
    if ~isKey(features,num2str(id))

        feat.timestamps=[];    
        feat.uvs_norm{cam_id}=[];     
        feat.uvs{cam_id}=[];    
        feat.featid=id;
        
        features(num2str(featid))=feat;
        
    end
    
    feat=features(num2str(featid));
        
    feat.timestamps=[feat.timestamps;timestamps];
    
    feat.uvs_norm{cam_id}=[feat.uvs_norm{cam_id};[u_n,v_n]];
    
    feat.uvs{cam_id}=[feat.uvs{cam_id};[u,v]];
        
    features(num2str(featid))=feat;
    
end
