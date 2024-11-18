function  err=computeError0(p1,p2,E)
    
    
    %E=V_2_Skew(T)*R;
    

    v1=[p1(1);p1(2);1];
    v2=[p2(1);p2(2);1];
    
    d1=v1'*E*v2;
    
    v1_=E*v2;
    
    v2_=E'*v1;
    
    d2=d1;
    
    s1=1/(v1_(1)^2+v1_(2)^2);
    
    s2=1/(v2_(1)^2+v2_(2)^2);
    
    
    err=max(d1*d1*s1,d2*d2*s2);
    
    err=sqrt(err);
   
end





function Skew=V_2_Skew(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end