function dRa_dtheta=D_Ra_D_theta(theta,a)
    
% the derivative of R(theta)*a with respect to theta 
%
%R=angleAxisToRotationMatrix(theta);
%
%Skew=Skew_symmetric(a);
%
%Jr_theta=J_r(theta);
    
dRa_dtheta=-angleAxisToRotationMatrix(theta)*Skew_symmetric(a)*J_r(theta);
     
    
end

