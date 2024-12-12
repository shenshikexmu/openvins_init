function  q=getInitQuaternionfromAcc_zfront(acc)


%有问题

norm_acc=norm([acc(1),acc(2),acc(3)]);


r31=acc(1)/norm_acc;

r32=acc(2)/norm_acc;

r33=acc(3)/norm_acc;


if r33==1
        
    R=[1,0,0;0,1,0;0,0,1];

elseif r33==-1

    R=[1,0,0;0,-1,0;0,0,-1];

else
    
    r23=0;
    
    r13=sqrt(1-r33^2);
    
    r21=-r31*r32/r22;
   
    r23=-r32*r33/r22;
    
    
    a=r33/r22;
    
    b=-r31/r22;
    
    if a*(r22*r33-r23*r32)+b*(r21*r32-r22*r31)>0
        
        r11=a;
        r13=b;
    else
        r11=-a;
        r13=-b;
    end
    
    
    R=[r11,r12,r13;...
       r21,r22,r23;...
       r31,r32,r33];
    
    
end


q=rotMat2qRichard(R)';



theta=acos(R(3,2))/pi*180;




end
