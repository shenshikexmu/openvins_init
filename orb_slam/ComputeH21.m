function H21 = ComputeH21(vP1, vP2)
    % Number of points
    N = length(vP1);

    % Initialize matrix A (2N x 9)
    A = zeros(2 * N, 9);

    for i = 1:N
        u1 = vP1(i, 1);  % x-coordinate of point in image 1
        v1 = vP1(i, 2);  % y-coordinate of point in image 1
        u2 = vP2(i, 1);  % x-coordinate of point in image 2
        v2 = vP2(i, 2);  % y-coordinate of point in image 2

        % Fill rows 2i-1 (indexing in MATLAB starts from 1)
        A(2*i-1, 1) = 0.0;
        A(2*i-1, 2) = 0.0;
        A(2*i-1, 3) = 0.0;
        A(2*i-1, 4) = -u1;
        A(2*i-1, 5) = -v1;
        A(2*i-1, 6) = -1;
        A(2*i-1, 7) = v2 * u1;
        A(2*i-1, 8) = v2 * v1;
        A(2*i-1, 9) = v2;

        % Fill rows 2i (indexing in MATLAB starts from 1)
        A(2*i, 1) = u1;
        A(2*i, 2) = v1;
        A(2*i, 3) = 1;
        A(2*i, 4) = 0.0;
        A(2*i, 5) = 0.0;
        A(2*i, 6) = 0.0;
        A(2*i, 7) = -u2 * u1;
        A(2*i, 8) = -u2 * v1;
        A(2*i, 9) = -u2;
    end

    % Perform SVD decomposition
    [~, ~, V] = svd(A);

    % Homography matrix is the last row of V (reshape it into a 3x3 matrix)
    H21 = reshape(V(:, end), 3, 3)';
end