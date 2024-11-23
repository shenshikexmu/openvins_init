function angleAxis=quaternionToAngleAxis(q)


sin_square_theta=q(2)*q(2)+q(3)*q(3)+q(4)*q(4);

if sin_square_theta>1e-9
    
    sin_theta=sqrt(sin_square_theta);
    
    cos_theta=q(1);
    
    if cos_theta<0
        
        two_theta=2*atan2(-sin_theta,-cos_theta);
    else
        two_theta=2*atan2(sin_theta,cos_theta);
    end
    k=two_theta/sin_theta;
    
    angleAxis(1,1)=q(2)*k;
    angleAxis(2,1)=q(3)*k;
    angleAxis(3,1)=q(4)*k;
else
    
    angleAxis(1,1)=q(2)*2;
    angleAxis(2,1)=q(3)*2;
    angleAxis(3,1)=q(4)*2;
    
end



end