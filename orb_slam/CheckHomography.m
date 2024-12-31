function [score, vbMatchesInliers] = CheckHomography(H21, H12, mvMatches12, mvKeys1, mvKeys2, sigma)

    % Number of matches
    N = size(mvMatches12, 1);

    % Extract homography matrix H21 and H12 elements
    h11 = H21(1, 1); h12 = H21(1, 2); h13 = H21(1, 3);
    h21 = H21(2, 1); h22 = H21(2, 2); h23 = H21(2, 3);
    h31 = H21(3, 1); h32 = H21(3, 2); h33 = H21(3, 3);

    h11inv = H12(1, 1); h12inv = H12(1, 2); h13inv = H12(1, 3);
    h21inv = H12(2, 1); h22inv = H12(2, 2); h23inv = H12(2, 3);
    h31inv = H12(3, 1); h32inv = H12(3, 2); h33inv = H12(3, 3);

    vbMatchesInliers = false(N, 1);  % Vector to store inliers
    score = 0;  % Initialize score

    % Threshold for inliers
    th = 5.991;

    % Precompute the inverse of sigma squared
    invSigmaSquare = 1.0 / (sigma * sigma);

    for i = 1:N
        bIn = true;

        % Extract keypoints for the pair of matches
        kp1 = mvKeys1(mvMatches12(i, 1), :);
        kp2 = mvKeys2(mvMatches12(i, 1), :);

        u1 = kp1(1); v1 = kp1(2);
        u2 = kp2(1); v2 = kp2(2);

        % Reprojection error in first image: x2in1 = H12 * x2
        w2in1inv = 1.0 / (h31inv * u2 + h32inv * v2 + h33inv);
        u2in1 = (h11inv * u2 + h12inv * v2 + h13inv) * w2in1inv;
        v2in1 = (h21inv * u2 + h22inv * v2 + h23inv) * w2in1inv;

        squareDist1 = (u1 - u2in1)^2 + (v1 - v2in1)^2;
        chiSquare1 = squareDist1 * invSigmaSquare;

        if chiSquare1 > th
            bIn = false;
        else
            score = score + (th - chiSquare1);
        end

        % Reprojection error in second image: x1in2 = H21 * x1
        w1in2inv = 1.0 / (h31 * u1 + h32 * v1 + h33);
        u1in2 = (h11 * u1 + h12 * v1 + h13) * w1in2inv;
        v1in2 = (h21 * u1 + h22 * v1 + h23) * w1in2inv;

        squareDist2 = (u2 - u1in2)^2 + (v2 - v1in2)^2;
        chiSquare2 = squareDist2 * invSigmaSquare;

        if chiSquare2 > th
            bIn = false;
        else
            score = score + (th - chiSquare2);
        end

        % Update inliers vector
        if bIn
            vbMatchesInliers(i) = true;
        else
            vbMatchesInliers(i) = false;
        end
    end
end
