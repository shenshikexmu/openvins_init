function outputImage = cv_resize(inputImage, dsize,fx,fy, interpolation)
% imageResize 缩放图像
%
% 输入:
%   inputImage - 输入图像（灰度或彩色）
%   dsize=[u, v]
%   fx x轴缩放因子
%   fy y轴缩放因子
%   interpolation 插值方式  INTER_NEAREST: 最近邻插值   INTER_LINEAR: 双线性插值（默认）
% 输出:
%   outputImage - 缩放后的图像

% 获取输入图像的尺寸
[rows, cols, channels] = size(inputImage);

if (dsize(1)~=0 && dsize(2)~=0)

    fy=dsize(1)/rows;
    fx=dsize(2)/cols;

end

dataType = class(inputImage);  % 确认类型

% 计算输出图像的尺寸
newRows = round(rows * fy);
newCols = round(cols * fx);

% 初始化输出图像
%outputImage = zeros(newRows, newCols, channels);   
outputImage = cast(zeros(newRows, newCols, channels),dataType);

% 定义插值系数
xCoeffs = linspace(0, 1, cols+1)';
yCoeffs = linspace(0, 1, rows+1);

% 遍历输出图像的每个像素
for i = 1:newRows
    for j = 1:newCols
        % 计算对应输入图像中的位置
        x = (j-0.5) / fx + 0.5;
        y = (i-0.5) / fy + 0.5;
        
        % 根据插值方法选择像素
        if strcmp(interpolation, 'INTER_NEAREST')
            % 最邻近插值
            x_nearest = round(x);
            y_nearest = round(y);
            
            % 确保索引在有效范围内
            x_nearest = max(min(x_nearest, cols), 1);
            y_nearest = max(min(y_nearest, rows), 1);
            
            % 为每个通道赋值
            for c = 1:channels
                outputImage(i, j, c) = inputImage(y_nearest, x_nearest, c);
            end

         elseif strcmp(interpolation, 'INTER_LINEAR')
        
             % 找到最近的四个输入图像像素
             x1 = floor(x);
             x2 = min(ceil(x), cols);
             y1 = floor(y);
             y2 = min(ceil(y), rows);
             
             % 确保索引在有效范围内
             if x1 < 1, x1 = 1; end
             if x2 < 1, x2 = 1; end
             if y1 < 1, y1 = 1; end
             if y2 < 1, y2 = 1; end
             
             % 计算双线性插值权重
             wx1 = 1 - (x - x1);
             wx2 = 1 - wx1;
             wy1 = 1 - (y - y1);
             wy2 = 1 - wy1;
         
             % 对每个通道进行插值
             for c = 1:channels
                 % 计算插值后的像素值
                 outputImage(i, j, c) = ...
                     cast(wy1 * (wx1 * inputImage(y1, x1, c) + wx2 * inputImage(y1, x2, c)) + ...
                           wy2 * (wx1 * inputImage(y2, x1, c) + wx2 * inputImage(y2, x2, c)),dataType);
             end
        else
            error('Unsupported interpolation method. Use "INTER_NEAREST" or "INTER_LINEAR".');
        end
    end
end


end