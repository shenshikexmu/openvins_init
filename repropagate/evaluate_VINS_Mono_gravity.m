function [residuals,Ji,Jj]=evaluate_VINS_Mono_gravity(Pi,Qi,Vi,Bai,Bgi,Pj,Qj,Vj,Baj,Bgj,imuPropagate,gravity)
   
 G=gravity;
    
delta_p=imuPropagate{1};
delta_q=imuPropagate{2};
delta_v=imuPropagate{3};
Ba=imuPropagate{4};
Bg=imuPropagate{5};
sum_dt=imuPropagate{6};
jacobian=imuPropagate{7};
covariance=imuPropagate{8};
precision=inv(covariance);

%[Q,R]=qr(precision);

dp_dba=jacobian(1:3,10:12);
dp_dbg=jacobian(1:3,13:15);
    
dq_dbg=jacobian(4:6,13:15);    
    
dv_dba=jacobian(7:9,10:12);
dv_dbg=jacobian(7:9,13:15);

dba=Bai-Ba;

dbg=Bgi-Bg;

 
corrected_delta_q=quaternProd(delta_q,  utility_deltaQ_VINS_Mono(dq_dbg*dbg));

corrected_delta_q=corrected_delta_q/norm(corrected_delta_q);

corrected_delta_v = delta_v + dv_dba * dba + dv_dbg * dbg;
    
corrected_delta_p = delta_p + dp_dba * dba + dp_dbg * dbg;  

residuals(1:3,1)=quatern2rotMat(Qi)'*(0.5*G*sum_dt*sum_dt+Pj-Pi-Vi*sum_dt)-corrected_delta_p;

residuals(4:6,1)= quaternionToAngleAxis( quaternProd( invQuaternion(corrected_delta_q), quaternProd(invQuaternion(Qi),Qj)) );
    
residuals(7:9,1)= quatern2rotMat(Qi)'*(G*sum_dt+Vj-Vi)-corrected_delta_v;

residuals(10:12,1)=Baj-Bai;

residuals(13:15,1)=Bgj-Bgi;


dqi_dthetai=-LeftMultiply( quaternProd(invQuaternion(Qj),Qi))*RightMultiply(corrected_delta_q);



J0=[-quatern2rotMat(Qi)',Skew_symmetric(quatern2rotMat(Qi)'*(0.5*G*sum_dt*sum_dt+Pj-Pi-Vi*sum_dt));...
    zeros(3,3),dqi_dthetai(2:4,2:4);...
    zeros(3,3),Skew_symmetric(quatern2rotMat(Qi)'*(G*sum_dt+Vj-Vi));...
    zeros(6,6)];

dqi_dbg=-LeftMultiply( quaternProd( quaternProd(invQuaternion(Qj),Qi), corrected_delta_q )  );

J1=[-quatern2rotMat(Qi)'*sum_dt,-dp_dba,-dp_dbg;...
    zeros(3,3),zeros(3,3),dqi_dbg(2:4,2:4)*dq_dbg;...
    -quatern2rotMat(Qi)',-dv_dba,-dv_dbg;...
    zeros(3,3),-eye(3),zeros(3,3);...
    zeros(3,3),zeros(3,3),-eye(3)];


dqj_dthetaj=LeftMultiply( quaternProd(quaternProd( invQuaternion(corrected_delta_q),invQuaternion(Qi)) , Qj)  );


J2=[quatern2rotMat(Qi)',zeros(3,3);...
    zeros(3,3), dqj_dthetaj(2:4,2:4)  ;...
    zeros(9,6)] ;

J3=[zeros(6,9);...
    quatern2rotMat(Qi)',zeros(3,6);...
    zeros(3,3),eye(3),zeros(3,3);...
    zeros(3,3),zeros(3,3),eye(3)];


%residuals_=residuals;

precision_dba=1;

precision_dbg=1;


for i=1:15
    
    residuals(i,1)=residuals(i,1)*sqrt(precision(i,i));
    
    J0(i,:)=J0(i,:)*sqrt(precision(i,i));
    
    J1(i,:)=J1(i,:)*sqrt(precision(i,i));
    
    J2(i,:)=J2(i,:)*sqrt(precision(i,i));
    
    J3(i,:)=J3(i,:)*sqrt(precision(i,i));
    
    
end

% precision_dba=500;
% 
% precision_dbg=1000;
% 
% 
% residuals(16:18,1)=dba*precision_dba;
% 
% residuals(19:21,1)=dbg*precision_dbg;


J0=[J0;zeros(6,6)];
J1=[J1;zeros(3,3),eye(3)*precision_dba,zeros(3,3);zeros(3,3),zeros(3,3),eye(3)*precision_dbg];
J2=[J2;zeros(6,6)];
J3=[J3;zeros(6,9)];

Ji=[J0,J1];

Jj=[J2,J3];


    
end



function R=RightMultiply(q)
    
R=[q(1),-q(2),-q(3),-q(4);...
   q(2), q(1), q(4),-q(3);...
   q(3),-q(4), q(1), q(2);...
   q(4), q(3),-q(2), q(1)];    
    
    
end



function R=LeftMultiply(q)
    
R=[q(1),-q(2),-q(3),-q(4);...
   q(2), q(1),-q(4), q(3);...
   q(3), q(4), q(1),-q(2);...
   q(4),-q(3), q(2), q(1)];    
    
    
end


function dRa_dq=D_Ra_D_q(q,a)

w=q(1);
v=[q(2);q(3);q(4)];
dRa_dq=2*[w*a+cross(v,a),v'*a*eye(3)+v*a'-a*v'-w*Skew_symmetric(a)];

end
