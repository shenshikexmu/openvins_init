function [vNormalizedPoints, T] = Normalize(vKeys)
    N = length(vKeys);  % Number of keypoints
    
    % Initialize variables
    meanX = 0;
    meanY = 0;
    
    % Calculate mean of the keypoints
    for i = 1:N
        meanX = meanX + vKeys(i,1);
        meanY = meanY + vKeys(i,2);
    end
    
    meanX = meanX / N;
    meanY = meanY / N;
    
    % Initialize normalized points and mean deviation variables
    vNormalizedPoints = zeros(N,2);
    meanDevX = 0;
    meanDevY = 0;
    
    % Normalize the points and compute the mean deviations
    for i = 1:N
        vNormalizedPoints(i,:) = [vKeys(i,1) - meanX, vKeys(i,2) - meanY];
        
        meanDevX = meanDevX + abs(vNormalizedPoints(i,1));
        meanDevY = meanDevY + abs(vNormalizedPoints(i,2));
    end
    
    meanDevX = meanDevX / N;
    meanDevY = meanDevY / N;
    
    % Compute scaling factors
    sX = 1.0 / meanDevX;
    sY = 1.0 / meanDevY;
    
    % Apply scaling to normalized points
%     for i = 1:N
%         vNormalizedPoints{i} = vNormalizedPoints{i} * [sX, sY];
%     end
    vNormalizedPoints=[vNormalizedPoints(:,1)*sX,vNormalizedPoints(:,2)*sY];
    
    % Construct the normalization matrix T
%     T = eye(3);
%     T(1,1) = sX;
%     T(2,2) = sY;
%     T(1,3) = -meanX * sX;
%     T(2,3) = -meanY * sY;
    T=[sX, 0,-meanX * sX;...
        0,sY,-meanY * sY;...
        0, 0,        1];
end
