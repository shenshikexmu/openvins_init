function dRa_dq=D_Ra_D_q(q,a)

w=q(1);
v=[q(2);q(3);q(4)];
dRa_dq=2*[w*a+cross(v,a),v'*a*eye(3)+v*a'-a*v'-w*Skew_symmetric(a)];

end