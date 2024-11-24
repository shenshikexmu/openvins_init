function q=utility_deltaQ_VINS_Mono(v)
    
    
q=[1;v(1)/2;v(2)/2;v(3)/2];
    
q=q/norm(q);
    
    
end