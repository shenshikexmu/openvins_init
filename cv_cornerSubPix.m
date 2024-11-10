function corners = cv_cornerSubPix(image, corners, winSize, zeroZone, criteria)
    % Ensure image is single channel (grayscale)
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    image = double(image);
    
    max_iters = 100;  % default max iterations
    eps = 0.01;       % default epsilon
    if criteria.type == "CV_TERMCRIT_ITER"
        max_iters = criteria.max_iter;
    elseif criteria.type == "CV_TERMCRIT_EPS"
        eps = criteria.epsilon;
    elseif criteria.type == "CV_TERMCRIT_ITER+EPS"
        eps = criteria.epsilon;
        max_iters = criteria.max_iter;
    end
    
    eps = max(eps, 0);
    eps = eps^2;  % use square of error in comparison operations
    max_iters = min(max(max_iters, 1), 100);  % clamp between 1 and 100

    % Initialize buffer for masks
    win_w = winSize(1) * 2 + 1;
    win_h = winSize(2) * 2 + 1;
    maskX = exp(-( (-winSize(1):winSize(1)) ).^2 / (winSize(1)^2));
    if winSize(1) == winSize(2)
        maskY = maskX;
    else
        maskY = exp(-( (-winSize(2):winSize(2)) ).^2 / (winSize(2)^2));
    end
    full_mask = maskX' * maskY;
    
    % Apply zero zone
    if zeroZone(1) >= 0 && zeroZone(2) >= 0
        for i = winSize(2) - zeroZone(2):winSize(2) + zeroZone(2)
            for j = winSize(1) - zeroZone(1):winSize(1) + zeroZone(1)
                full_mask(i + 1, j + 1) = 0;
            end
        end
    end
    
    % Iterate over each corner
    for pt_i = 1:size(corners, 1)
        cT = corners(pt_i, :);
        cI = cT;
        iter = 0;
        err = Inf;
        
        while iter < max_iters && err > eps
            % Get subpixel neighborhood and adjust mask size
            [sub_pix, actual_win_w, actual_win_h] = getSubPix(image, cI, win_w, win_h);

            
            [m_sub_pix,n_sub_pix]=size(sub_pix);
            if (m_sub_pix<2 || n_sub_pix<2)
                cI=cT;
                break;
            end
            
            % Resize mask to match actual window size
            mask = imresize(full_mask, [actual_win_h, actual_win_w], 'bilinear');

                       
            % Compute gradients
            [gx, gy] = gradient(sub_pix);
            
            % Compute Hessian components and vector
            a = 0; b = 0; c = 0;
            bb1 = 0; bb2 = 0;
            for i = 1:actual_win_h
                py = i - floor(actual_win_h / 2) - 1;
                for j = 1:actual_win_w
                    px = j - floor(actual_win_w / 2) - 1;
                    gxx = gx(i,j)^2 * mask(i,j);
                    gxy = gx(i,j) * gy(i,j) * mask(i,j);
                    gyy = gy(i,j)^2 * mask(i,j);
                    a = a + gxx;
                    b = b + gxy;
                    c = c + gyy;
                    bb1 = bb1 + gxx * px + gxy * py;
                    bb2 = bb2 + gxy * px + gyy * py;
                end
            end
            
            % Solve for dx, dy
            detA = a * c - b^2;
            if abs(detA) < 1e-10
                break;
            end
            dx = (c * bb1 - b * bb2) / detA;
            dy = (a * bb2 - b * bb1) / detA;
            
            cI2 = cI + [dx, dy];
            err = (cI2(1) - cI(1))^2 + (cI2(2) - cI(2))^2;
            cI = cI2;
            iter = iter + 1;
        end
        
        if norm(cI - cT) > max(winSize)
            cI = cT;  % Poor convergence, revert to initial point
        end
        
        corners(pt_i, :) = cI;  % Store result
    end
end

function [sub_pix, actual_win_w, actual_win_h] = getSubPix(image, center, win_w, win_h)
    % Get subpixel region from the image around center point
    % Ensure the coordinates stay within image bounds
    
    [img_h, img_w] = size(image);  % Get image dimensions
    x = round(center(1));
    y = round(center(2));
    
    % Compute window boundaries, ensuring they are within the image bounds
    x_min = max(1, x - floor(win_w / 2));
    x_max = min(img_w, x + floor(win_w / 2));
    y_min = max(1, y - floor(win_h / 2));
    y_max = min(img_h, y + floor(win_h / 2));
    
    % Extract the sub-pixel region
    sub_pix = image(y_min:y_max, x_min:x_max);
    
    % Get actual window size (may be smaller than original if close to borders)
    actual_win_w = x_max - x_min + 1;
    actual_win_h = y_max - y_min + 1;
end



% function corners = cornerSubPix(image, corners, winSize, zeroZone, criteria)
%     % Ensure image is single channel (grayscale)
%     if size(image, 3) == 3
%         image = rgb2gray(image);
%     end
%     image = double(image);
%     
%     max_iters = 100;  % default max iterations
%     eps = 0.01;       % default epsilon
%     if criteria.type == "CV_TERMCRIT_ITER"
%         max_iters = criteria.max_iter;
%     elseif criteria.type == "CV_TERMCRIT_EPS"
%         eps = criteria.epsilon;
%     elseif criteria.type == "CV_TERMCRIT_ITER+EPS"
%         eps = criteria.epsilon;
%         max_iters = criteria.max_iter;
%     end
%     
%     eps = max(eps, 0);
%     eps = eps^2;  % use square of error in comparison operations
%     max_iters = min(max(max_iters, 1), 100);  % clamp between 1 and 100
% 
%     % Initialize buffer for masks
%     win_w = winSize(1) * 2 + 1;
%     win_h = winSize(2) * 2 + 1;
%     maskX = exp(-( (-winSize(1):winSize(1)) ).^2 / (winSize(1)^2));
%     if winSize(1) == winSize(2)
%         maskY = maskX;
%     else
%         maskY = exp(-( (-winSize(2):winSize(2)) ).^2 / (winSize(2)^2));
%     end
%     mask = maskX' * maskY;
%     
%     % Apply zero zone
%     if zeroZone(1) >= 0 && zeroZone(2) >= 0
%         for i = winSize(2) - zeroZone(2):winSize(2) + zeroZone(2)
%             for j = winSize(1) - zeroZone(1):winSize(1) + zeroZone(1)
%                 mask(i + 1, j + 1) = 0;
%             end
%         end
%     end
%     
%     % Iterate over each corner
%     for pt_i = 1:size(corners, 1)
%         cT = corners(pt_i, :);
%         cI = cT;
%         iter = 0;
%         err = Inf;
%         
%         while iter < max_iters && err > eps
%             % Get subpixel neighborhood
%             sub_pix = getSubPix(image, cI, win_w, win_h);
%             
%             % Compute gradients
%             [gx, gy] = gradient(sub_pix);
%             
%             % Compute Hessian components and vector
%             a = 0; b = 0; c = 0;
%             bb1 = 0; bb2 = 0;
%             for i = 1:win_h
%                 py = i - winSize(2) - 1;
%                 for j = 1:win_w
%                     px = j - winSize(1) - 1;
%                     gxx = gx(i,j)^2 * mask(i,j);
%                     gxy = gx(i,j) * gy(i,j) * mask(i,j);
%                     gyy = gy(i,j)^2 * mask(i,j);
%                     a = a + gxx;
%                     b = b + gxy;
%                     c = c + gyy;
%                     bb1 = bb1 + gxx * px + gxy * py;
%                     bb2 = bb2 + gxy * px + gyy * py;
%                 end
%             end
%             
%             % Solve for dx, dy
%             detA = a * c - b^2;
%             if abs(detA) < 1e-10
%                 break;
%             end
%             dx = (c * bb1 - b * bb2) / detA;
%             dy = (a * bb2 - b * bb1) / detA;
%             
%             cI2 = cI + [dx, dy];
%             err = (cI2(1) - cI(1))^2 + (cI2(2) - cI(2))^2;
%             cI = cI2;
%             iter = iter + 1;
%         end
%         
%         if norm(cI - cT) > max(winSize)
%             cI = cT;  % Poor convergence, revert to initial point
%         end
%         
%         corners(pt_i, :) = cI;  % Store result
%     end
% end
% 
% function sub_pix = getSubPix(image, center, win_w, win_h)
%     % Get subpixel region from the image around center point
%     x = round(center(1));
%     y = round(center(2));
%     sub_pix = image(y - floor(win_h/2):y + floor(win_h/2), x - floor(win_w/2):x + floor(win_w/2));
% end
