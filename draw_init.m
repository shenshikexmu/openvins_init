function []=draw_init(features,R1,T1,R2,T2)


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


ids_all = keys(features);

for i = 1:length(ids_all)

    id = ids_all{i}; % 获取当前键

    feat = features(id);

    if ~isempty( feat.p_FinA)

        plot3(feat.p_FinA(1),feat.p_FinA(2),feat.p_FinA(3),'r*');

    end

end

hold off;


end