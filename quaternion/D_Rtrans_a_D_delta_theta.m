function dRtrans_a_ddelta_theta=D_Rtrans_a_D_delta_theta(q,a)

% the derivative of R(q)*a with respect to q 
% q(1) is Q.w
    
dRtrans_a_ddelta_theta=D_Rtrans_a_D_q(q,a)*D_q_D_delta_theta(q);

end