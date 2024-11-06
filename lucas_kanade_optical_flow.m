function [u, v] = lucas_kanade_optical_flow(I1, I2, window_size)
    % I1, I2: 两帧图像
    % window_size: 卷积窗口的大小，通常为 5 或 7
    
    % 转换为双精度以进行计算
    I1 = double(I1);
    I2 = double(I2);
    
    % 计算图像梯度
    [Ix, Iy] = gradient(I1);
    It = I2 - I1; % 时间梯度
    
    % 初始化光流 u, v
    [rows, cols] = size(I1);
    u = zeros(rows, cols);
    v = zeros(rows, cols);
    
    % 设置窗口半径
    half_window = floor(window_size / 2);
    
    % 遍历每个像素，计算光流
    for i = half_window+1:rows-half_window
        for j = half_window+1:cols-half_window
            % 提取当前窗口区域
            Ix_window = Ix(i-half_window:i+half_window, j-half_window:j+half_window);
            Iy_window = Iy(i-half_window:i+half_window, j-half_window:j+half_window);
            It_window = It(i-half_window:i+half_window, j-half_window:j+half_window);
            
            % 构造矩阵 A 和 B
            A = [sum(Ix_window(:).^2), sum(Ix_window(:).*Iy_window(:));
                 sum(Ix_window(:).*Iy_window(:)), sum(Iy_window(:).^2)];
            B = [-sum(Ix_window(:).*It_window(:));
                 -sum(Iy_window(:).*It_window(:))];
            
            % 计算光流 u 和 v
            if det(A) > 1e-5 % 确保 A 是非奇异的
                flow = A \ B;
                u(i, j) = flow(1);
                v(i, j) = flow(2);
            end
        end
    end
end
