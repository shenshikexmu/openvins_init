function R=Rodrigues_formula(w)
    
norm_w=norm(w);
    
if norm_w==0
    
    R=eye(3);

else
    
    Skew_w=Skew_symmetric(w);
    
    R=eye(3)+sin(norm_w)/norm_w*Skew_w+(1-cos(norm_w))/norm_w^2*Skew_w*Skew_w;
 
end
  
end



function Skew=Skew_symmetric(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end