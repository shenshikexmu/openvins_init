function Q=angleAxis2Quaternion(angleAxis)


R=angleAxisToRotationMatrix(angleAxis);

Q=rotMat2qRichard_nolo(R);


end


function R=angleAxisToRotationMatrix(angleAxis)


theta2=angleAxis(1)*angleAxis(1)+angleAxis(2)*angleAxis(2)+angleAxis(3)*angleAxis(3);



if theta2>1e-9
    
    theta=sqrt(theta2);
    costheta=cos(theta);
    sintheta=sin(theta);
    theta_inverse=1/theta;
    n=[angleAxis(1);angleAxis(2);angleAxis(3)]*theta_inverse;
    
    
    R=costheta*eye(3)+(1-costheta)*n*n'+sintheta*V_2_Skew(n);
    
else
    
    R=eye(3);
    
end



end


function Skew=V_2_Skew(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end
