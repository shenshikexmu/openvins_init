function Jr= rightJacobian_of_SO3(w)
    
    
norm_w=norm(w);    
    
    
if norm_w==0
    
    Jr=eye(3);

else
    
    Skew_w=V_2_Skew(w);
    
    Jr=eye(3)-(1-cos(norm_w))/norm_w^2*Skew_w+(norm_w-sin(norm_w))/norm_w^3*Skew_w*Skew_w;

end
    
    
end




function Skew=V_2_Skew(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end