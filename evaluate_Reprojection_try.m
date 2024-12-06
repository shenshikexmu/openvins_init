function [E,H_dz_dG_I_p,H_dz_dG_I_angleAxis,H_dz_dG_p_f_k]=evaluate_Reprojection_try(G_I_p,G_I_angleAxis,G_p_f_k,camK,camD,camR,camT,pts_temp)

G_I_R=angleAxisToRotationMatrix(G_I_angleAxis);

G_c_R=G_I_R*camR;

G_c_T=G_I_p+G_I_R*camT;

C_p_f=G_c_R'*(G_p_f_k-G_c_T);


[uv_dist,H_dz_dzn]= distort_cv(C_p_f(1:2)/C_p_f(3), camK,camD);


E=(uv_dist-pts_temp)';


d_uv_norm_d_C_p_f=[1/C_p_f(3),0,-C_p_f(1)/C_p_f(3)^2;...
                   0,1/C_p_f(3),-C_p_f(2)/C_p_f(3)^2];
                   
d_C_p_f_d_G_p_f_k=G_c_R';

d_C_p_f_d_G_I_p=-G_c_R';


A=G_p_f_k-G_c_T;

d_C_p_f_d_G_I_angleAxis=camR'*D_Rtrans_a_D_theta(G_I_angleAxis,A);

H_dz_dG_p_f_k=H_dz_dzn*d_uv_norm_d_C_p_f*d_C_p_f_d_G_p_f_k;

H_dz_dG_I_p=H_dz_dzn*d_uv_norm_d_C_p_f*d_C_p_f_d_G_I_p;


H_dz_dG_I_angleAxis=H_dz_dzn*d_uv_norm_d_C_p_f*d_C_p_f_d_G_I_angleAxis;

end