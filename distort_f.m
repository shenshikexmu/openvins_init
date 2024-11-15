function uv_dist = distort_f(uv_norm, camK,camD)
    
%    camera_values =[fx, fy, cx, cy, k1, k2, p1, p2]
%    camK=[fx, 0,cx;...
%           0,fy,cy;...
%           0, 0, 1];
%    camD=[k1,k2,p1,p2];

    % 计算径向畸变
    r_2 = uv_norm(1)^2 + uv_norm(2)^2;
    r_4 = r_2^2;
    
    % 计算畸变后的坐标
    x1 = uv_norm(1) * (1 + camD(1) * r_2 + camD(2) * r_4) + 2 * camD(3) * uv_norm(1) * uv_norm(2) + camD(4) * (r_2 + 2 * uv_norm(1)^2);
    y1 = uv_norm(2) * (1 + camD(1) * r_2 + camD(2) * r_4) + cam_d(3) * (r_2 + 2 * uv_norm(2)^2) + 2 * camD(4) * uv_norm(1) * uv_norm(2);

    % 返回畸变后的点
    uv_dist = [camK(1,1) * x1 + camK(1,3), camK(2,2) * y1 + camK(2,3)];
end