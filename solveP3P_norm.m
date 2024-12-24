function [Rall,Tall]=solveP3P_norm(object_points, image_points_norm)

Rall=[];
Tall=[];
    
n=0;  

for i=1:size(object_points)
   
    n=n+1;
    PP(n,:)=[image_points_norm(i,1),image_points_norm(i,2),object_points(i,1),object_points(i,2),object_points(i,3)];
    
    if n==3
        break
    end
end

Vx=[PP(1,1),PP(1,2),1];
Vy=[PP(2,1),PP(2,2),1];
Vz=[PP(3,1),PP(3,2),1];
Pa=PP(1,3:5);
Pb=PP(2,3:5);
Pc=PP(3,3:5);

cos_yz=  Vy*Vz'/(norm(Vy)*norm(Vz)) ;
cos_xz=  Vx*Vz'/(norm(Vx)*norm(Vz)) ;
cos_xy=  Vx*Vy'/(norm(Vx)*norm(Vy)) ;

p=2*cos_yz;
q=2*cos_xz;
r=2*cos_xy;
    
AB2=(Pa-Pb)*(Pa-Pb)';
AC2=(Pa-Pc)*(Pa-Pc)';
BC2=(Pb-Pc)*(Pb-Pc)';

a=BC2/AB2;
b=AC2/AB2;


% D=[ -a  , 1-a ,  0 , -p , a*r , 1;...
%     1-b , -b  , -q ,  0 , b*r , 1 ];


a0= -2*b+ b^2+ a^2+ 1- b*r^2*a+ 2*b*a- 2*a;

a1= -2*b*q*a- 2*a^2*q+ b*r^2*q*a- 2*q+ 2*b*q+ 4*a*q+ p*b*r+ b*r*p*a- b^2*r*p;

a2= q^2+ b^2*r^2- b*p^2- q*p*b*r+ b^2*p^2- b*r^2*a+ 2- 2*b^2- a*b*r*p*q+ 2*a^2- 4*a- 2*q^2*a+ q^2*a^2;

a3=-b^2*r*p+ b*r*p*a- 2*a^2*q+ q*p^2*b+ 2*b*q*a+ 4*a*q+ p*b*r- 2*b*q- 2*q;

a4= 1- 2*a+ 2*b+ b^2- b*p^2+ a^2- 2*b*a;

%solve('a0*x^4+a1*x^3+a2*x^2+a3*x+a4=0',x);

% syms x
% eqn=a0*x^4+a1*x^3+a2*x^2+a3*x+a4==0;
% vpasolve(eqn)

%x=sovle('a0*x^4+a1*x^3+a2*x^2+a3*x+a4=0','x')

xv=roots([a0,a1,a2,a3,a4]);

n=0;

for i=1:4
    
    if imag(xv(i,1)) ==0
        
        n=n+1;
        xV(n,1)=real(xv(i,1));
           
    end
     
end



for i=1:size(xV,1)
    
    x=xV(i,1);

    b0=b*(p^2*a- p^2+ b*p^2+ p*q*r- q*a*r*p+ a*r^2- r^2- b*r^2)^2;

    b1=((1-a-b)*x^2+ (q*a-q)*x+ 1-a+b)*(...
        ...
        (a^2*r^3+ 2*b*r^3*a- b*r^5*a- 2*a*r^3+ r^3+ b^2*r^3- 2*r^3*b )*x^3+ ...
        ...
        (p*r^2+ p*a^2*r^2- 2*b*r^3*q*a+ 2*r^3*b*q- 2*r^3*q- 2*p*a*r^2- 2*p*r^2*b+ ...
         r^4*p*b+ 4*a*r^3*q+ b*q*a*r^5- 2*r^3*a^2*q+ 2*r^2*p*b*a+ b^2*r^2*p- r^4*p*b^2 )*x^2+ ...
        ...
        (r^3*q^2+ r^5*b^2+ r*p^2*b^2- 4*a*r^3- 2*a*r^3*q^2+ r^3*q^2*a^2+...
         2*a^2*r^3- 2*b^2*r^3- 2*p^2*b*r+ 4*p*a*r^2*q+ 2*a*p^2*r*b-...
         2*a*r^2*q*b*p- 2*p^2*a*r+ r*p^2- b*r^5*a+ 2*p*r^2*b*q+ ...
         r*p^2*a^2- 2*p*q*r^2+ 2*r^3- 2*r^2*p*a^2*q- r^4*q*b*p )*x+ ...
        ...
         4*a*r^3*q+ p*r^2*q^2+ 2*p^3*b*a- 4*p*a*r^2- ...
         2*r^3*b*q- 2*p^2*q*r- 2*b^2*r^2*p+ r^4*p*b+ 2*p*a^2*r^2- ...
         2*r^3*a^2*q- 2*p^3*a+ p^3*a^2+ 2*p*r^2+ p^3+ 2*b*r^3*q*a+ ...
         2*q*p^2*b*r+4*q*a*r*p^2-2*p*a*r^2*q^2-2*p^2*a^2*r*q+ ...
         p*a^2*r^2*q^2- 2*r^3*q- 2*p^3*b+ p^3*b^2- 2*p^2*b*r*q*a ...
         );
     
     y=b1/b0;
     
     yV(i,1)=y;


end


% for i=1:size(xV,1)                           %
%                                              
%     x=xV(i);
%     y=yV(i);
% 
%     p1=(1-a)*y^2-a*x^2-p*y+a*r*x*y +1;
% 
%     p2=(1-b)*x^2-b*y^2-q*x+b*r*x*y +1;
%     
%     PPPP(1:2,i)=[p1;p2];
%     
% end

n=0;


for i=1:size(xV,1)
    
    x=xV(i,1);
    y=yV(i,1);
    
    lz=sqrt((AC2/(1+x^2-x*q)+BC2/(1+y^2-y*p))/2);
    
    %Z22=bc2/(1+y^2-y*p);
    
    lx=x*lz;
    
    ly=y*lz;
    
    
    PA=lx*Vx'/norm(Vx);
    
    PB=ly*Vy'/norm(Vy);
    
    PC=lz*Vz'/norm(Vz);
    
    
    %[[norm(PA-PB), norm(PA-PC) ,norm(PB-PC) ]' ,   [norm(Pa-Pb) ,norm(Pa-Pc) ,norm(Pb-Pc)]']
     
    %[PA,PB,PC]

    
    PABCmean=(PA+PB+PC)/3;
    
    Pa=PP(1,3:5)';
    Pb=PP(2,3:5)';
    Pc=PP(3,3:5)';

    Pabcmean=(Pa+Pb+Pc)/3;
    
    
    PA_=PA-PABCmean;
    PB_=PB-PABCmean;
    PC_=PC-PABCmean;
    
    
    Pa_=Pa-Pabcmean;
    Pb_=Pb-Pabcmean;
    Pc_=Pc-Pabcmean;
    
    
    W=PA_*Pa_'+PB_*Pb_'+PC_*Pc_';
    
    
    [UU,SS,VV]=svd(W);
    
    if det(UU)<-0.9
        UU=[UU(:,1),UU(:,2),-UU(:,3)];
    end
        
    if det(VV)<-0.9
        VV=[VV(:,1),VV(:,2),-VV(:,3)];
    end    
        
    
    %[W,UU,SS,VV]
    
    
    R=UU*VV';
    
    
%     PABC_=[PA_,PB_,PC_];
%     
%     Pabc_=[Pa_,Pb_,Pc_];
 
    T=(PA-R*Pa+PB-R*Pb+PC-R*Pc)/3;
        
    %RTP3P=[R,T]
    
    n=n+1;
    
    Rall{n}=R;
    
    Tall{n}=T;
    
    
end

    
    
    
    
end


