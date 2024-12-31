function x3D = Triangulate(kp1, kp2, P1, P2)
    % Initialize matrix A (4x4) for the triangulation
    A = zeros(4, 4);
    
    % Set up the rows of matrix A
    A(1, :) = kp1(1) * P1(3, :) - P1(1, :);  % First row
    A(2, :) = kp1(2) * P1(3, :) - P1(2, :);  % Second row
    A(3, :) = kp2(1) * P2(3, :) - P2(1, :);  % Third row
    A(4, :) = kp2(2) * P2(3, :) - P2(2, :);  % Fourth row
    
    % Perform Singular Value Decomposition (SVD) on matrix A
    [~, ~, V] = svd(A);
    
    % Extract the last row of V (which corresponds to the smallest singular value)
    x3D = V(:, end);
    
    % Normalize the 3D point
    x3D = x3D(1:3) / x3D(4);
end
