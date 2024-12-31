function [score,vbMatchesInliers] = CheckFundamental(F21, mvMatches12, mvKeys1, mvKeys2, sigma)
    % Number of matches
    N = size(mvMatches12, 1);

    % Extract fundamental matrix F21 elements
    f11 = F21(1, 1); f12 = F21(1, 2); f13 = F21(1, 3);
    f21 = F21(2, 1); f22 = F21(2, 2); f23 = F21(2, 3);
    f31 = F21(3, 1); f32 = F21(3, 2); f33 = F21(3, 3);

    vbMatchesInliers = false(N, 1);  % Initialize inliers vector
    score = 0;  % Initialize score

    % Threshold for inliers
    th = 3.841;
    thScore = 5.991;

    % Precompute the inverse of sigma squared
    invSigmaSquare = 1.0 / (sigma * sigma);

    for i = 1:N
        bIn = true;

        % Extract keypoints for the pair of matches
        kp1 = mvKeys1(mvMatches12(i, 1), :);
        kp2 = mvKeys2(mvMatches12(i, 1), :);

        u1 = kp1(1); v1 = kp1(2);
        u2 = kp2(1); v2 = kp2(2);

        % Reprojection error in the second image: l2 = F21 * x1
        a2 = f11 * u1 + f12 * v1 + f13;
        b2 = f21 * u1 + f22 * v1 + f23;
        c2 = f31 * u1 + f32 * v1 + f33;

        num2 = a2 * u2 + b2 * v2 + c2;

        squareDist1 = num2^2 / (a2^2 + b2^2);

        chiSquare1 = squareDist1 * invSigmaSquare;

        if chiSquare1 > th
            bIn = false;
        else
            score = score + (thScore - chiSquare1);
        end

        % Reprojection error in the first image: l1 = x2' * F21
        a1 = f11 * u2 + f21 * v2 + f31;
        b1 = f12 * u2 + f22 * v2 + f32;
        c1 = f13 * u2 + f23 * v2 + f33;

        num1 = a1 * u1 + b1 * v1 + c1;

        squareDist2 = num1^2 / (a1^2 + b1^2);

        chiSquare2 = squareDist2 * invSigmaSquare;

        if chiSquare2 > th
            bIn = false;
        else
            score = score + (thScore - chiSquare2);
        end

        % Update inliers vector
        if bIn
            vbMatchesInliers(i) = true;
        else
            vbMatchesInliers(i) = false;
        end
    end
end
