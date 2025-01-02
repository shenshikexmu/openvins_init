function [vbMatchesInliers, score, H21] = FindHomography(mvMatches12, mvKeys1, mvKeys2, mvSets, mMaxIterations, mSigma)
% Number of putative matches
N = size(mvMatches12, 1);

% Normalize coordinates
[vPn1, T1] = Normalize(mvKeys1);
[vPn2, T2] = Normalize(mvKeys2);
T2inv = inv(T2);

% Best Results variables
score = 0.0;
vbMatchesInliers = false(N, 1);

% Iteration variables
vPn1i = zeros(8, 2);
vPn2i = zeros(8, 2);
H21i = [];
H12i = [];
vbCurrentInliers = false(N, 1);
currentScore = 0;

% Perform all RANSAC iterations and save the solution with the highest score
for it = 1:mMaxIterations
    % Select a minimum set
%     for j = 1:8
%         idx = mvSets{it}(j);
% 
%         vPn1i(j, :) = vPn1(mvMatches12{idx}(1), :);
%         vPn2i(j, :) = vPn2(mvMatches12{idx}(2), :);
%     end
    vPn1i=vPn1(mvMatches12(mvSets(:,it)), :);
    vPn2i=vPn2(mvMatches12(mvSets(:,it)), :);

    Hn = ComputeH21(vPn1i, vPn2i);
    H21i = T2inv * Hn * T1;
    H12i = inv(H21i);

    [currentScore, vbCurrentInliers] = CheckHomography(H21i, H12i, mvMatches12, mvKeys1, mvKeys2, mSigma);

    if currentScore > score
        H21 = H21i;
        vbMatchesInliers = vbCurrentInliers;
        score = currentScore;
    end
end

end
