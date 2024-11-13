function corners = cv_FAST(image, threshold, nonmaxSuppression)
    % 输入:
    % image - 输入的灰度图像，像素值为 0-255
    % threshold - 判断像素与中心像素差异的阈值
    % nonmaxSuppression - 是否启用非极大值抑制 (true/false)
    
    % 输出:
    % corners - 角点坐标及响应值的矩阵 [x, y, response]
    
    % 获取图像大小
    [rows, cols] = size(image);
    
    % 定义圆形邻域的 16 个像素相对位置（顺时针方向）
    circleOffsets = [0, 3; 1, 3; 2, 2; 3, 1; 3, 0; 3, -1; 2, -2; 1, -3; 0, -3; -1, -3; -2, -2; -3, -1; -3, 0; -3, 1; -2, 2; -1, 3];
    
    % 初始化角点响应矩阵
    responseMatrix = zeros(rows, cols);
    
    % 遍历每个像素点，排除边界
    for r = 4:rows-3
        for c = 4:cols-3
            % 获取中心像素
            
            p = image(r, c);
            
            % 提取邻域 16 个像素
            circlePixels = zeros(16, 1);
            for i = 1:16
                circlePixels(i) = image(r + circleOffsets(i, 1), c + circleOffsets(i, 2));
            end
            
            % 对比像素值并统计比中心像素明显更亮或更暗的像素数量
            brighter = circlePixels > p + threshold;
            darker = circlePixels < p - threshold;
            
            % 如果有连续 9 个像素更亮或更暗，认定为角点
            if any(conv([brighter; brighter(1:8)], ones(9, 1), 'valid') == 9) || ...
               any(conv([darker; darker(1:8)], ones(9, 1), 'valid') == 9)
                % 计算响应值（这里简单定义为亮和暗像素差异的最大值）
                responseMatrix(r, c) = max(abs(circlePixels - double(p)));
            end
        end
    end
    
    % 非极大值抑制（去掉局部邻域内较弱的角点）
    if nonmaxSuppression
        responseMatrix = nonmax_suppression(responseMatrix);
    end
    
    % 获取角点的坐标和响应值
    [cornerRows, cornerCols] = find(responseMatrix > 0);
    cornerResponses = responseMatrix(responseMatrix > 0);
    
    % 将角点坐标和响应值组合成输出矩阵
    corners = [cornerCols, cornerRows, cornerResponses];

    corners = sortrows(corners, 2);
end

function suppressed = nonmax_suppression(responseMatrix)
    % 非极大值抑制，保留局部最大值，抑制其他较弱的响应
    suppressed = responseMatrix;
    [rows, cols] = size(responseMatrix);
    
    for r = 2:rows-1
        for c = 2:cols-1
            localPatch = responseMatrix(r-1:r+1, c-1:c+1);
            if responseMatrix(r, c) < max(localPatch(:))
                suppressed(r, c) = 0; % 抑制非极大值
            end
        end
    end
end
