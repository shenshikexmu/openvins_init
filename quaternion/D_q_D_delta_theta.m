function  dq_ddelta_theta=D_q_D_delta_theta(q)
    
% the derivative of Q with respect to delta_theta 
% q(1) is Q.w
    
dq_ddelta_theta=0.5*[-q(2),-q(3),-q(4);...
                      q(1),-q(4), q(3);...
                      q(4), q(1),-q(2);...
                     -q(3), q(2), q(1)];
    
end
