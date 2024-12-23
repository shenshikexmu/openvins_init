function dRtans_a_dtheta=D_Rtrans_a_D_theta(theta,a)
    
% the derivative of Rtans(theta)*a with respect to theta 
%
%R=angleAxisToRotationMatrix(theta);
%
%Skew=Skew_symmetric(a);
%
%Jr_theta=J_r(theta);
    
dRtans_a_dtheta=Skew_symmetric(angleAxisToRotationMatrix(theta)'*a)*J_r(theta);
     
    
end

