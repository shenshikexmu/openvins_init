function [uv_dist,H_dz_dzn] = distort_cv(uv_norm, camK,camD)

    [uv_dist,H_dz_dzn] = distort_f(uv_norm, camK,camD);


end



function [uv_dist,H_dz_dzn] = distort_f(uv_norm, camK,camD)
    
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
    y1 = uv_norm(2) * (1 + camD(1) * r_2 + camD(2) * r_4) + camD(3) * (r_2 + 2 * uv_norm(2)^2) + 2 * camD(4) * uv_norm(1) * uv_norm(2);

    % 返回畸变后的点
    uv_dist = [camK(1,1) * x1 + camK(1,3), camK(2,2) * y1 + camK(2,3)];
    
    [H_dz_dzn, H_dz_dzeta] = compute_distort_jacobian(uv_norm, camK,camD);
    
end



function [H_dz_dzn, H_dz_dzeta] = compute_distort_jacobian(uv_norm, camK,camD)
    
    % Calculate distorted coordinates for radial distortion
    r_2 = uv_norm(1)^2 + uv_norm(2)^2;
    r_4 = r_2^2;

    x = uv_norm(1);
    y = uv_norm(2);
    x_2 = x^2;
    y_2 = y^2;
    x_y = x * y;

    % Jacobian of distorted pixel to normalized pixel (H_dz_dzn)
    H_dz_dzn = zeros(2, 2);
    H_dz_dzn(1, 1) = camK(1,1) * ((1 + camD(1) * r_2 + camD(2) * r_4) + ...
                                 (2 * camD(1) * x_2 + 4 * camD(2) * x_2 * r_2) + ...
                                 2 * camD(3) * y + (2 * camD(4) * x + 4 * camD(4) * x));
    H_dz_dzn(1, 2) = camK(1,1) * (2 * camD(1) * x_y + 4 * camD(2) * x_y * r_2 + ...
                                 2 * camD(3) * x + 2 * camD(4) * y);
    H_dz_dzn(2, 1) = camK(2,2) * (2 * camD(1) * x_y + 4 * camD(2) * x_y * r_2 + ...
                                 2 * camD(3) * x + 2 * camD(4) * y);
    H_dz_dzn(2, 2) = camK(2,2) * ((1 + camD(1) * r_2 + camD(2) * r_4) + ...
                                 (2 * camD(1) * y_2 + 4 * camD(2) * y_2 * r_2) + ...
                                 2 * camD(4) * x + (2 * camD(3) * y + 4 * camD(3) * y));

    % Calculate distorted coordinates for radtan model
    x1 = x * (1 + camD(1) * r_2 + camD(2) * r_4) + 2 * camD(3) * x * y + ...
         camD(4) * (r_2 + 2 * x_2);
    y1 = y * (1 + camD(1) * r_2 + camD(2) * r_4) + ...
         camD(3) * (r_2 + 2 * y_2) + 2 * camD(4) * x * y;

    % Compute the Jacobian in respect to the intrinsics (H_dz_dzeta)
    H_dz_dzeta = zeros(2, 8);
    H_dz_dzeta(1, 1) = x1;
    H_dz_dzeta(1, 3) = 1;
    H_dz_dzeta(1, 5) = camK(1,1) * x * r_2;
    H_dz_dzeta(1, 6) = camK(1,1) * x * r_4;
    H_dz_dzeta(1, 7) = 2 * camK(1,1) * x * y;
    H_dz_dzeta(1, 8) = camK(1,1) * (r_2 + 2 * x_2);
    H_dz_dzeta(2, 2) = y1;
    H_dz_dzeta(2, 4) = 1;
    H_dz_dzeta(2, 5) = camK(2,2) * y * r_2;
    H_dz_dzeta(2, 6) = camK(2,2) * y * r_4;
    H_dz_dzeta(2, 7) = camK(2,2) * (r_2 + 2 * y_2);
    H_dz_dzeta(2, 8) = 2 * camK(2,2) * x * y;
end

