function [nGood, vP3D, vbGood, parallax,allError] = CheckRT(R, t, vKeys1, vKeys2, vMatches12, vbMatchesInliers, K, th2)
    
% Calibration parameters
fx = K(1, 1);
fy = K(2, 2);
cx = K(1, 3);
cy = K(2, 3);

N = length(vKeys1);
vbGood = false(N, 1);
vP3D = zeros(3, N);  % 3D points
vCosParallax = [];   % To store the cos of parallax
allError=0;

% Camera 1 Projection Matrix K[I|0]
P1 = [K, zeros(3, 1)];
O1 = [0; 0; 0]; % Camera center for the first camera

% Camera 2 Projection Matrix K[R|t]
P2 = K * [R, t];
O2 = -R' * t; % Camera center for the second camera

nGood = 0; % Counter for good points

for i = 1:length(vMatches12)
    if ~vbMatchesInliers(i)
        continue;
    end
    
    % Get the keypoints from the match
    kp1 = vKeys1(vMatches12(i, 1),:);
    kp2 = vKeys2(vMatches12(i, 1),:);
    
    % Triangulate the point
    p3dC1 = Triangulate(kp1, kp2, P1, P2); % Returns 3D point in camera 1 frame
    
    if any(~isfinite(p3dC1)) % If any coordinate is NaN or Inf, skip it
        vbGood(vMatches12(i, 1)) = false;
        continue;
    end

    % Compute parallax
    normal1 = p3dC1 - O1;
    dist1 = norm(normal1);
    normal2 = p3dC1 - O2;
    dist2 = norm(normal2);
    
    cosParallax = dot(normal1, normal2) / (dist1 * dist2);

    % Check depth in front of the first camera (ensure valid point)
    if p3dC1(3) <= 0 && cosParallax < 0.99998
        continue;
    end

    % Check depth in front of the second camera (ensure valid point)
    p3dC2 = R * p3dC1 + t;
    if p3dC2(3) <= 0 && cosParallax < 0.99998
        continue;
    end

    % Reprojection error in the first image
    invZ1 = 1 / p3dC1(3);
    im1x = fx * p3dC1(1) * invZ1 + cx;
    im1y = fy * p3dC1(2) * invZ1 + cy;
    squareError1 = (im1x - kp1(1))^2 + (im1y - kp1(2))^2;

%     if squareError1 > th2
%         continue;
%     end

    % Reprojection error in the second image
    invZ2 = 1 / p3dC2(3);
    im2x = fx * p3dC2(1) * invZ2 + cx;
    im2y = fy * p3dC2(2) * invZ2 + cy;
    squareError2 = (im2x - kp2(1))^2 + (im2y - kp2(2))^2;

%     if squareError2 > th2
%         continue;
%     end
    allError=allError+squareError1+squareError2;

    % Store the 3D point and parallax
    vCosParallax = [vCosParallax; cosParallax];
    vP3D(:,vMatches12(i, 1)) = p3dC1;

    if  p3dC2(3)>0 &&  p3dC1(3)>0
        nGood = nGood + 1;
    end

    % If cosParallax is low, consider it a good match
    if cosParallax < 0.99998
        vbGood(vMatches12(i, 1)) = true;
    end
end

% If there are good points, compute the parallax
if nGood > 0
    vCosParallax = sort(vCosParallax);
    idx = min(20, length(vCosParallax));
    parallax = acos(vCosParallax(idx)) * 180 / pi;
else
    parallax = 0;
end

end
