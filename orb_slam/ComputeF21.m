function F21 = ComputeF21(vP1, vP2)
    % Number of points
    N = size(vP1, 1);

    % Initialize matrix A (N x 9)
    A = zeros(N, 9);

    for i = 1:N
        u1 = vP1(i, 1);  % x-coordinate of point in image 1
        v1 = vP1(i, 2);  % y-coordinate of point in image 1
        u2 = vP2(i, 1);  % x-coordinate of point in image 2
        v2 = vP2(i, 2);  % y-coordinate of point in image 2

        % Fill matrix A
        A(i, 1) = u2 * u1;
        A(i, 2) = u2 * v1;
        A(i, 3) = u2;
        A(i, 4) = v2 * u1;
        A(i, 5) = v2 * v1;
        A(i, 6) = v2;
        A(i, 7) = u1;
        A(i, 8) = v1;
        A(i, 9) = 1;
    end

    % Perform SVD decomposition
    [~, ~, V] = svd(A);

    % Fundamental matrix Fpre is the last row of V (reshape it into a 3x3 matrix)
    Fpre = reshape(V(:, end), 3, 3)';

    % Perform another SVD on Fpre to enforce rank-2 constraint
    [Uf, Sf, Vf] = svd(Fpre);

    % Set the smallest singular value to zero to enforce rank-2 constraint
    Sf(3, 3) = 0;

    % Reconstruct the fundamental matrix with rank-2 constraint
    F21 = Uf * Sf * Vf';
end
