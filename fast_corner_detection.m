function keypoints = fast_corner_detection(img, threshold)
    % FAST corner detection algorithm in MATLAB
    % img: Grayscale input image
    % threshold: Threshold for corner detection

    % Convert to grayscale if not already
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    
    % Define offsets for Bresenham circle (16 pixels)
    circle_offsets = [
        0, -3;
        1, -3;
        2, -2;
        3, -1;
        3, 0;
        3, 1;
        2, 2;
        1, 3;
        0, 3;
        -1, 3;
        -2, 2;
        -3, 1;
        -3, 0;
        -3, -1;
        -2, -2;
        -1, -3;
    ];

    % Initialize the keypoints array
    keypoints = [];

    % Get image dimensions
    [rows, cols] = size(img);
    
    % Loop over each pixel, excluding borders
    for y = 4:rows-3
        for x = 4:cols-3
            % Get the intensity of the center pixel
            center_pixel = img(y, x);
            
            % Get the pixel intensities on the circle
            circle_pixels = zeros(16, 1);
            for i = 1:16
                circle_pixels(i) = img(y + circle_offsets(i, 2), x + circle_offsets(i, 1));
            end
            
            % Compare pixel intensities with threshold
            brighter = circle_pixels > (center_pixel + threshold);
            darker = circle_pixels < (center_pixel - threshold);
            
            % Check if there are at least 9 consecutive brighter or darker pixels
            if any(conv([brighter; brighter], ones(9, 1), 'valid') == 9) || ...
               any(conv([darker; darker], ones(9, 1), 'valid') == 9)
                % Mark this as a keypoint
                keypoints = [keypoints; x, y]; %#ok<AGROW>
            end
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end
