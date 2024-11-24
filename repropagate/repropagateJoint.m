function imuPropagate_1_3=repropagateJoint(imuPropagate_1_2,imuPropagate_2_3)

%function imuPropagate_1_3=repropagateJoint(imuPropagate_1_2,imuPropagate_2_3,imuPropagate_1_3_)


P_1_2=imuPropagate_1_2{1};
Q_1_2=imuPropagate_1_2{2};
V_1_2=imuPropagate_1_2{3};
t_1_2=imuPropagate_1_2{6};
% Jac_1_2=imuPropagate_1_2{7};
% Cov_1_2=imuPropagate_1_2{8};

P_2_3=imuPropagate_2_3{1};
Q_2_3=imuPropagate_2_3{2};
V_2_3=imuPropagate_2_3{3};
t_2_3=imuPropagate_2_3{6};
% Jac_2_3=imuPropagate_2_3{7};
% Cov_2_3=imuPropagate_2_3{8};


% P_1_3_=imuPropagate_1_3_{1};
% Q_1_3_=imuPropagate_1_3_{2};
% V_1_3_=imuPropagate_1_3_{3};
% Jac_1_3_=imuPropagate_1_3_{7};
% Cov_1_3_=imuPropagate_1_3_{8};



R_1_2=quatern2rotMat(Q_1_2);


P_1_3=P_1_2+R_1_2*P_2_3+V_1_2*t_2_3;

Q_1_3=quaternProd(Q_1_2,Q_2_3);

V_1_3=V_1_2+R_1_2*V_2_3;


imuPropagate_1_3{1}=P_1_2+R_1_2*P_2_3+V_1_2*t_2_3;
imuPropagate_1_3{2}=quaternProd(Q_1_2,Q_2_3);
imuPropagate_1_3{3}=V_1_2+R_1_2*V_2_3;
imuPropagate_1_3{4}=imuPropagate_1_2{4};
imuPropagate_1_3{5}=imuPropagate_1_2{5};
imuPropagate_1_3{6}=imuPropagate_1_2{6}+imuPropagate_2_3{6};


% Jac_1_3=Jac_2_3*Jac_1_2;  jac还有问题




end