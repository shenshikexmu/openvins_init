function dRtrans_a_dq=D_Rtrans_a_D_q(q,a)

% the derivative of R(q)*a with respect to q 
% q(1) is Q.w
    
w=q(1);
v=[q(2);q(3);q(4)];
dRtrans_a_dq=2*[w*a+cross(a,v),v'*a*eye(3)+v*a'-a*v'+w*Skew_symmetric(a)];

end