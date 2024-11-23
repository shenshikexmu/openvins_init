function uv_norm=cv_undistortPoints(uv_dist, camK,camD)

%    camera_values =[fx, fy, cx, cy, k1, k2, p1, p2]
%    camK=[fx, 0,cx;...
%           0,fy,cy;...
%           0, 0, 1];
%    camD=[k1,k2,p1,p2];

xd=[(uv_dist(1)-camK(1,3))/camK(1,1);(uv_dist(2)-camK(2,3))/camK(2,2)];

k1 = camD(1);
k2 = camD(2);
k3=0;
p1 = camD(3);
p2 = camD(4);

x = xd; 				% initial guess

for kk=1:20

    r_2 = sum(x.^2);
    k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
    delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
    p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
    x = (xd - delta_x)./(ones(2,1)*k_radial);

end
    
uv_norm=[x(1),x(2)];    
    
end





% function [x] = comp_distortion_oulu(xd,k)
% 
%     
% k1 = k(1);
% k2 = k(2);
% k3 = k(5);
% p1 = k(3);
% p2 = k(4);
% 
% x = xd; 				% initial guess
% 
% for kk=1:20
% 
%     r_2 = sum(x.^2);
%     k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
%     delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
%     p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
%     x = (xd - delta_x)./(ones(2,1)*k_radial);
% 
% end
%     
% end

















