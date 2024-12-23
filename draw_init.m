function []=draw_init(features,map_camera_times,R1,T1,R2,T2,cam_id1,cam_id2,frame1,frame2)



timestamp1=map_camera_times(frame1,1)-map_camera_times(frame1,3);
timestamp2=map_camera_times(frame2,1)-map_camera_times(frame2,3);



figure

R{1}=R1;
T{1}=T1;
R{2}=R2;
T{2}=T2;

colors = jet(size(R,2));

scale=0.5;

img_x_len=0.2*scale;
img_y_len=0.127*scale;

for i=1:size(R,2)

    G_camR_k=R{i};  
    G_camT_k=T{i};
    
    p1=G_camT_k-G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
    p2=G_camT_k-G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
    p3=G_camT_k+G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
    p4=G_camT_k+G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
    
    lines=[p1,p2,p3,p4,p1];
    
    plot3(lines(1,:),lines(2,:),lines(3,:),'-', 'Color', colors(i,:));
    hold on;
    if i==1

        text(G_camT_k(1),G_camT_k(2),G_camT_k(3),'0');
        hold on

        plot3([0,1]*0.1*scale,[0,0],[0,0],'r-',...
              [0,0],[0,1]*0.1*scale,[0,0],'g-',...
              [0,0],[0,0],[0,1]*0.1*scale,'b-' );

    end
     
end

axis equal;

Len=10;


ids_all = keys(features);

for i = 1:length(ids_all)

    id = ids_all{i}; % 获取当前键

    feat = features(id);

    if ~isempty( feat.p_FinA)

        plot3(feat.p_FinA(1),feat.p_FinA(2),feat.p_FinA(3),'r*');


        flag1=0; 

        if size(feat.timestamps{cam_id1},1)>1
    
            for j=1:size(feat.timestamps{cam_id1},1)
    
                if timestamp1==feat.timestamps{cam_id1}(j)

                    flag1=1;
                    
                    uv_n=feat.uvs_norm{cam_id1}(j,:);
                    V1=R1*[uv_n(1);uv_n(2);1];

                    

    
                end
    
            end
    
        end
    
        flag2=0;

        if size(feat.timestamps{cam_id2},1)>1
    
            for j=1:size(feat.timestamps{cam_id2},1)
    
                if timestamp2==feat.timestamps{cam_id2}(j)

                    flag2=1;


                    uv_n=feat.uvs_norm{cam_id2}(j,:);
                    V2=R2*[uv_n(1);uv_n(2);1];

                  

                    
                   
                end
    
            end
    
        end

        if flag1==1 && flag2==1

            hold on

            plot3([T1(1),T1(1)+V1(1)*Len],[T1(2),T1(2)+V1(2)*Len],[T1(3),T1(3)+V1(3)*Len],'r-',...
                  [T2(1),T2(1)+V2(1)*Len],[T2(2),T2(2)+V2(2)*Len],[T2(3),T2(3)+V2(3)*Len],'g-');

        end

    end

end

hold off;


end