function [R1, R2, t] = DecomposeE(E)
    % Perform Singular Value Decomposition of the essential matrix
    [U, ~, V] = svd(E);

    % Extract translation vector
    t = U(:, 3); % The third column of U is the translation vector
    t = t / norm(t); % Normalize the translation vector

    % Define the W matrix
    W = [  0, -1,  0;
            1,  0,  0;
            0,  0,  1];

    % Compute the first rotation matrix R1
    R1 = U * W * V';
    if det(R1) < 0
        R1 = -R1; % Ensure that the determinant is positive
    end

    % Compute the second rotation matrix R2
    R2 = U * W' * V';
    if det(R2) < 0
        R2 = -R2; % Ensure that the determinant is positive
    end

end