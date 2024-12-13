function R_GtoI = gram_schmidt(gravity_inI)
    
    
    % This will find an orthogonal vector to gravity which is our local z-axis
    % We need to ensure we normalize after each one such that we obtain unit vectors
    z_axis = gravity_inI / norm(gravity_inI);
    x_axis = [];
    y_axis = [];
    e_1 = [1.0; 0.0; 0.0];
    e_2 = [0.0; 1.0; 0.0];
    inner1 = dot(e_1, z_axis) / norm(z_axis); % norm(z_axis) is 1 since it's already normalized
    inner2 = dot(e_2, z_axis) / norm(z_axis); % same here
    
    if abs(inner1) < abs(inner2)
        x_axis = cross(z_axis, e_1);
        x_axis = x_axis / norm(x_axis);
        y_axis = cross(z_axis, x_axis);
        y_axis = y_axis / norm(y_axis);
    else
        x_axis = cross(z_axis, e_2);
        x_axis = x_axis / norm(x_axis);
        y_axis = cross(z_axis, x_axis);
        y_axis = y_axis / norm(y_axis);
    end

    % Rotation from our global (where gravity is only along the z-axis) to the local one
    R_GtoI = [x_axis, y_axis, z_axis];
end