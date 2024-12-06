function Jr_theta=J_r(theta)
    
 norm_theta=norm(theta);
    
 Jr_theta=eye(3)- (1-cos(norm_theta))/norm_theta^2*Skew_symmetric(theta)+(norm_theta-sin(norm_theta))/norm_theta^3*Skew_symmetric(theta)*Skew_symmetric(theta);
    
end