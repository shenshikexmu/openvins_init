function [u, v] = lucas_kanade_optical_flow_points(I1, I2, points, window_size)
    % 输入:
    % I1: 第一帧图像 (灰度)
    % I2: 第二帧图像 (灰度)
    % points: 特征点列表 (Nx2的矩阵，每行是一个[x, y]点的坐标)
    % window_size: 窗口大小
    
    % 输出:
    % u: 特征点的水平光流分量 (N x 1 向量)
    % v: 特征点的垂直光流分量 (N x 1 向量)
    
    % 计算图像梯度
    [Ix, Iy] = gradient(double(I1)); % X方向和Y方向的梯度
    It = double(I2) - double(I1);    % 时间上的梯度 (帧间差异)
    
    num_points = size(points, 1);
    u = zeros(num_points, 1);
    v = zeros(num_points, 1);
    
    half_w = floor(window_size / 2);
    
    % 对每个特征点计算光流
    for k = 1:num_points
        % 获取当前特征点的坐标
        x = round(points(k, 1));
        y = round(points(k, 2));
        
        % 检查特征点是否在图像范围内
        if x-half_w < 1 || x+half_w > size(I1, 2) || y-half_w < 1 || y+half_w > size(I1, 1)
            % 如果超出边界，跳过这个点
            continue;
        end
        
        % 提取窗口内的梯度
        Ix_window = Ix(y-half_w:y+half_w, x-half_w:x+half_w);
        Iy_window = Iy(y-half_w:y+half_w, x-half_w:x+half_w);
        It_window = It(y-half_w:y+half_w, x-half_w:x+half_w);
        
        % 将梯度展开为向量形式
        Ix_vec = Ix_window(:);
        Iy_vec = Iy_window(:);
        It_vec = -It_window(:); % 注意这里负号
        
        % 构建A矩阵和b向量 (A*u = b)
        A = [Ix_vec, Iy_vec];
        b = It_vec;
        
        % 检查A是否可逆
        if rank(A) == 2
            % 使用最小二乘法求解
            nu = pinv(A' * A) * A' * b;
            u(k) = nu(1);
            v(k) = nu(2);
        else
            % 如果A不可逆，光流设为0
            u(k) = 0;
            v(k) = 0;
        end
    end
end
