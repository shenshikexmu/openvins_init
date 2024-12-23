function invJr= invRightJacobian_of_SO3(w)
    
    
norm_w=norm(w);    
    
    
if norm_w==0
    
    invJr=eye(3);

else
    
    Skew_w=Skew_symmetric(w);
    
    invJr=eye(3)+1/2*Skew_w+(1/norm_w^2+(1+cos(norm_w)/(2*norm_w*sin(norm_w))))* Skew_w* Skew_w;
    

end
    
    
end




function Skew=Skew_symmetric(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end