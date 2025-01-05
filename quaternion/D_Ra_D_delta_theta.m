function dRa_ddelta_theta=D_Ra_D_delta_theta(q,a)

% the derivative of R(q)*a with respect to q 
% q(1) is Q.w
    
dRa_ddelta_theta=D_Ra_D_q(q,a)*D_q_D_delta_theta(q);

end