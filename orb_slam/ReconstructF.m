function [success, R21, t21,maxGood,vP3D, vbTriangulated] = ReconstructF(vbMatchesInliers, F21, K, mvMatches12, mvKeys1, mvKeys2,minParallax, minTriangulated,mSigma2)

success = false;

% Number of matches
N = sum(vbMatchesInliers);

% Compute Essential Matrix from Fundamental Matrix
E21 = K' * F21 * K;

% Recover the 4 motion hypotheses (rotation and translation)
[R1, R2, t] = DecomposeE(E21);

% Set two possible values for translation
t1 = t;
t2 = -t;

% Initialize variables for 3D reconstruction and triangulation status
%     vP3D1 = []; vP3D2 = []; vP3D3 = []; vP3D4 = [];
%     vbTriangulated1 = []; vbTriangulated2 = []; vbTriangulated3 = []; vbTriangulated4 = [];
%     parallax1 = 0; parallax2 = 0; parallax3 = 0; parallax4 = 0;

% Check RT hypotheses
[nGood1, vP3D1, vbTriangulated1, parallax1] = CheckRT(R1, t1, mvKeys1, mvKeys2, mvMatches12, vbMatchesInliers, K, 4.0 * mSigma2);
[nGood2, vP3D2, vbTriangulated2, parallax2] = CheckRT(R2, t1, mvKeys1, mvKeys2, mvMatches12, vbMatchesInliers, K, 4.0 * mSigma2);
[nGood3, vP3D3, vbTriangulated3, parallax3] = CheckRT(R1, t2, mvKeys1, mvKeys2, mvMatches12, vbMatchesInliers, K, 4.0 * mSigma2);
[nGood4, vP3D4, vbTriangulated4, parallax4] = CheckRT(R2, t2, mvKeys1, mvKeys2, mvMatches12, vbMatchesInliers, K, 4.0 * mSigma2);

% Find the best hypothesis based on the number of good points
maxGood = max([nGood1, nGood2, nGood3, nGood4]);

R21 = []; t21 = []; vP3D=[]; vbTriangulated=[]; % Initialize output variables

% Minimum number of good triangulated points
nMinGood = max(floor(0.9 * N), minTriangulated);

% Check for similar hypotheses
nsimilar = sum([nGood1, nGood2, nGood3, nGood4] > 0.7 * maxGood);

% Reject if not enough triangulated points or there are too many similar hypotheses
if maxGood < nMinGood || nsimilar > 1
    success = false;
    return;
end

% Select the best hypothesis based on parallax
if maxGood == nGood1
    if parallax1 > minParallax
        vP3D = vP3D1;
        vbTriangulated = vbTriangulated1;
        R21 = R1;
        t21 = t1;
        success = true;
        return;
    end
elseif maxGood == nGood2
    if parallax2 > minParallax
        vP3D = vP3D2;
        vbTriangulated = vbTriangulated2;
        R21 = R2;
        t21 = t1;
        success = true;
        return;
    end
elseif maxGood == nGood3
    if parallax3 > minParallax
        vP3D = vP3D3;
        vbTriangulated = vbTriangulated3;
        R21 = R1;
        t21 = t2;
        success = true;
        return;
    end
elseif maxGood == nGood4
    if parallax4 > minParallax
        vP3D = vP3D4;
        vbTriangulated = vbTriangulated4;
        R21 = R2;
        t21 = t2;
        success = true;
        return;
    end
end


end
