function []=draw_init_camera(R,T,n)



colors = jet(10);

scale=2;

img_x_len=0.2*scale;
img_y_len=0.127*scale;



G_camR_k=R;  
G_camT_k=T;

p1=G_camT_k-G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
p2=G_camT_k-G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
p3=G_camT_k+G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
p4=G_camT_k+G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;

lines=[p1,p2,p3,p4,p1];

plot3(lines(1,:),lines(2,:),lines(3,:),'-', 'Color', colors(n,:));
hold on;



end
    
       

