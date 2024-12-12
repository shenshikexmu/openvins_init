function  q=getInitQuaternionfromAcc_xfront(acc)




norm_acc=norm([acc(1),acc(2),acc(3)]);


r31=acc(1)/norm_acc;

r32=acc(2)/norm_acc;

r33=acc(3)/norm_acc;


if r31==1
        
    R=[0,0,-1;0,1,0;1,0,0];

elseif r31==-1

    R=[0,0,1;0,0,1;-1,0,0];

else
    
    r21=0;
    
    r11=sqrt(1-r31^2);
    
    r12=-r31*r32/r11;
   
    r13=-r31*r33/r11;
    
    
    a=r33/r11;
    
    b=-r32/r11;
    
    if a*(r11*r33-r13*r31)+b*(r31*r12-r11*r32)>0
        
        r22=a;
        r23=b;
    else
        r22=-a;
        r23=-b;
    end
    
    
    R=[r11,r12,r13;...
       r21,r22,r23;...
       r31,r32,r33];
    
    
end


q=rotMat2qRichard(R);



theta=acos(R(3,2))/pi*180;




end
