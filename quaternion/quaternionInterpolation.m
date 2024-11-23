function q=quaternionInterpolation(q1,q2,a,b)
    


invq1= invQuaternion(q1);
    
q_delta=quaternProd(invq1,q2);
    
angleAxis_delta=quaternionToAngleAxis(q_delta);    


%Q_delta_backN=angleAxis2Quaternion(angleAxis_delta/backN);

q_delta_inter = axisAngle2quatern2(angleAxis_delta*b);

q=quaternProd(q1,q_delta_inter);
    
    
    
end