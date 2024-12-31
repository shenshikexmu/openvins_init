function [vbMatchesInliers, score, F21] = FindFundamental(mvMatches12, mvKeys1, mvKeys2, mvSets, mMaxIterations, mSigma)
% Number of putative matches
N = size(mvMatches12, 1);

% Normalize coordinates
[vPn1, T1] = Normalize(mvKeys1);
[vPn2, T2] = Normalize(mvKeys2);
T2t = T2';

% Best Results variables
score = 0.0;
vbMatchesInliers = false(N,1);

% Iteration variables
vPn1i = zeros(8, 2);
vPn2i = zeros(8, 2);
F21i = [];
vbCurrentInliers = false(N,1);
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
    vPn1i=vPn1(mvMatches12(mvSets(:,1)), :);
    vPn2i=vPn2(mvMatches12(mvSets(:,1)), :);

    % Compute Fundamental matrix
    Fn = ComputeF21(vPn1i, vPn2i);

    % Compute F21 using normalization transformation
    F21i = T2t * Fn * T1;

    % Check the score and update if necessary
    [currentScore,vbCurrentInliers] = CheckFundamental(F21i, mvMatches12, mvKeys1, mvKeys2, mSigma);

    if currentScore > score
        F21 = F21i;
        vbMatchesInliers = vbCurrentInliers;
        score = currentScore;
    end
end

end