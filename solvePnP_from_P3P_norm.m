function [R,T]=solvePnP_from_P3P_norm(object_points, image_points_norm)


n=0;

for i=1:size(object_points)
   
    n=n+1;
    PP(n,:)=[image_points_norm(i,1),image_points_norm(i,2),object_points(i,1),object_points(i,2),object_points(i,3)];

end

if n>3

    [R3,T3]=solveP3P_norm(object_points, image_points_norm);

end

min=inf;

for i=1:size(R3,2)
    
    error=0;
    
    for j=1:size(PP,1)
        
         UVtemp=[R3{i},T3{i}]*[PP(j,3:5)';1];
         
         error=error+norm(UVtemp(1:2)/UVtemp(3)-PP(j,1:2)');
        
    end
    
    if error < min 
        
        min=error;
        
        R_=R3{i};
        
        T_=T3{i};
        
    end
 
end


[R_,T_]=My_GaussNewton_PnP(PP,R_,T_);


R=R_';

T=-R_'*T_;



end



function [R,T]=My_GaussNewton_PnP(PP,Rinit,Tinit)

q=rotMat2qRichard(Rinit);

phi=quaternionToAngleAxis(q);


X=[phi;Tinit];


for i=1:size(PP,1)
     
     UV(i*2-1,1)=PP(i,1);
        
     UV(i*2-0,1)=PP(i,2);
    
end


for i=1:6
    
    [UVn,Jacbi]=f_PnP(X,PP);
    
    %UV-UVn
    
    X=X+inv(Jacbi'*Jacbi)*Jacbi'*(UV-UVn);
    
end

R=angleAxisToRotationMatrix(X(1:3));

T=X(4:6);

end



function [UVn,Jacbi]=f_PnP(X,PP)

angle=norm(X(1:3));

if angle==0
    
    angle_axis=[1;0;0];
    
    Jl=eye(3);
    
else
    
    angle_axis=X(1:3)/angle;
    
    Jl=sin(angle)/angle*eye(3)+(1-sin(angle)/angle)*angle_axis*angle_axis'+(1-cos(angle))/angle*Skew_symmetric(angle_axis); 
    
end
    
R=angleAxisToRotationMatrix(X(1:3));

T=X(4:6);

UVn=zeros(size(PP,1)*2,1);

Jacbi=zeros(size(PP,1)*2,6);

for i=1:size(PP,1)
    
    P=PP(i,3:5)';
    
    Vtemp=R*P+T;
    
    UVn(i*2-1,1)= Vtemp(1)/ Vtemp(3);
    
    UVn(i*2-0,1)= Vtemp(2)/ Vtemp(3);
    
    d_UV_d_Vtemp=[1/Vtemp(3),     0    ,- Vtemp(1)/ Vtemp(3)^2;...
                       0    ,1/Vtemp(3),- Vtemp(2)/ Vtemp(3)^2];
                   
    d_Vtemp_d_X=[ -Skew_symmetric(R*P)*Jl,eye(3)];
    
    Jacbi(i*2-1:i*2,:)=d_UV_d_Vtemp*d_Vtemp_d_X;
    
end





end




























