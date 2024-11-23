function q = axisAngle2quatern2(angleAxis)
   
    
    angle=norm(angleAxis);
    if angle==0
        
        q = [1 0 0 0];
    
    
    else
  
        axis=angleAxis/angle;

        q0 = cos(angle/2);
        q1 = axis(1)*sin(angle/2);
        q2 = axis(2)*sin(angle/2);
        q3 = axis(3)*sin(angle/2); 
        q = [q0 q1 q2 q3];
    end
end