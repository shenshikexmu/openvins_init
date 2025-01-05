function dRa_dq=D_Ra_D_q(q,a)

% the derivative of R(q)*a with respect to q 
% q(1) is Q.w
    
w=q(1);
v=[q(2);q(3);q(4)];
dRa_dq=2*[w*a+cross(v,a),v'*a*eye(3)+v*a'-a*v'-w*Skew_symmetric(a)];

end