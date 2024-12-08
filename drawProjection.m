function []=drawProjection(x_I_k,G_p_f,camR,camT)


figure

for i=1:size(G_p_f,2)

   G_p_f_k=G_p_f(:,i);
   plot3(G_p_f_k(1),G_p_f_k(2),G_p_f_k(3),'r*');
   hold on;

end


colors = jet(size(x_I_k,2));

scale=0.5;

img_x_len=0.2*scale;
img_y_len=0.127*scale;

for i=1:size(x_I_k,2)

    
    x_I_k_=x_I_k(:,i);
    
    G_p_k=x_I_k_(5:7);
    
    G_q_k=x_I_k_(1:4);
    
    G_R_k=quatern2rotMat(G_q_k);
    
    G_camR_k=G_R_k*camR;
    
    G_camT_k=G_p_k+G_R_k*camT;
    
    p1=G_camT_k-G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
    p2=G_camT_k-G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
    p3=G_camT_k+G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
    p4=G_camT_k+G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
    
    lines=[p1,p2,p3,p4,p1];
    
    plot3(lines(1,:),lines(2,:),lines(3,:),'-', 'Color', colors(i,:));
    hold on;
    if i==1
        text(G_camT_k(1),G_camT_k(2),G_camT_k(3),'0');

    end
     
end


axis equal;

hold off;




end