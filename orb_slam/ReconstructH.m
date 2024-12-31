function [success,R21, t21,vP3D, vbTriangulated] = ReconstructH(vbMatchesInliers, H21, K,mvMatches12, mvKeys1, mvKeys2, minParallax, minTriangulated,mSigma2)
    % Number of matches
    N = sum(vbMatchesInliers);
    
    % Compute the inverse of K
    invK = inv(K);
    
    % Compute A = inv(K) * H21 * K
    A = invK * H21 * K;

    % SVD decomposition of A
    [U, w, Vt] = svd(A, 'econ');
    V = Vt';

    % Check the determinant of U and Vt to ensure the validity of the decomposition
    s = det(U) * det(Vt);

    d1 = w(1,1);
    d2 = w(2,2);
    d3 = w(3,3);

    % If the ratio between d1, d2, and d3 is too small, return false
    if (d1 / d2 < 1.00001) || (d2 / d3 < 1.00001)
        success = false;
        return;
    end

    % Prepare vectors for storing rotation matrices, translation vectors, and normal vectors
    vR = cell(8, 1);
    vt = cell(8, 1);
    vn = cell(8, 1);

    % Compute values for rotation and translation hypotheses based on the Faugeras method
    aux1 = sqrt((d1^2 - d2^2) / (d1^2 - d3^2));
    aux3 = sqrt((d2^2 - d3^2) / (d1^2 - d3^2));
    x1 = [aux1, aux1, -aux1, -aux1];
    x3 = [aux3, -aux3, aux3, -aux3];

    aux_stheta = sqrt((d1^2 - d2^2) * (d2^2 - d3^2)) / ((d1 + d3) * d2);
    ctheta = (d2^2 + d1 * d3) / ((d1 + d3) * d2);
    stheta = [aux_stheta, -aux_stheta, -aux_stheta, aux_stheta];

    for i = 1:4
        Rp = eye(3);
        Rp(1,1) = ctheta;
        Rp(1,3) = -stheta(i);
        Rp(3,1) = stheta(i);
        Rp(3,3) = ctheta;

        R = s * U * Rp * Vt;
        vR{i} = R;

        tp = [x1(i); 0; -x3(i)] * (d1 - d3);
        t = U * tp;
        vt{i} = t / norm(t);

        np = [x1(i); 0; x3(i)];
        n = V * np;
        if n(3) < 0
            n = -n;
        end
        vn{i} = n;
    end

    % Case for d' = -d2
    aux_sphi = sqrt((d1^2 - d2^2) * (d2^2 - d3^2)) / ((d1 - d3) * d2);
    cphi = (d1 * d3 - d2^2) / ((d1 - d3) * d2);
    sphi = [aux_sphi, -aux_sphi, -aux_sphi, aux_sphi];

    for i = 1:4
        Rp = eye(3);
        Rp(1,1) = cphi;
        Rp(1,3) = sphi(i);
        Rp(2,2) = -1;
        Rp(3,1) = sphi(i);
        Rp(3,3) = -cphi;

        R = s * U * Rp * Vt;
        vR{4 + i} = R;

        tp = [x1(i); 0; x3(i)] * (d1 + d3);
        t = U * tp;
        vt{4 + i} = t / norm(t);

        np = [x1(i); 0; x3(i)];
        n = V * np;
        if n(3) < 0
            n = -n;
        end
        vn{4 + i} = n;
    end

    % Initialize variables for the best solution
    bestGood = 0;
    secondBestGood = 0;
    bestSolutionIdx = -1;
    bestParallax = -1;
    bestP3D = [];
    bestTriangulated = [];

    % Check each of the 8 hypotheses
    for i = 1:8
        parallaxi = 0;
        vP3Di = [];
        vbTriangulatedi = [];
        [nGood, vP3Di, vbTriangulatedi, parallaxi] = CheckRT(vR{i}, vt{i}, mvKeys1, mvKeys2, mvMatches12, vbMatchesInliers, K, 4.0 * mSigma2);


        if nGood > bestGood
            secondBestGood = bestGood;
            bestGood = nGood;
            bestSolutionIdx = i;
            bestParallax = parallaxi;
            bestP3D = vP3Di;
            bestTriangulated = vbTriangulatedi;
        elseif nGood > secondBestGood
            secondBestGood = nGood;
        end
    end

    % If the best solution has sufficient parallax and triangulated points, accept it
    if secondBestGood < 0.75 * bestGood && bestParallax >= minParallax && bestGood > minTriangulated && bestGood > 0.9 * N
        R21 = vR{bestSolutionIdx};
        t21 = vt{bestSolutionIdx};
        vP3D = bestP3D;
        vbTriangulated = bestTriangulated;
        success = true;
    else
        success = false;
        R21 = [];
        t21 = [];
        vP3D = [];
        vbTriangulated =[];
    end
end
