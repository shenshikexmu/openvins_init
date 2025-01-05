function [R,T]=optimaization_Fundamental_R_T(pts1_n, pts2_n,Rinit,Tinit,vP3D)


a=zeros(6+size(pts1_n,1),1);


q=rotMat2qRichard(Rinit);

phi=quaternionToAngleAxis(q);


a0(1:6,1)=[phi;Tinit];

for i=1:size(pts1_n,1)

    a0(6+i*3-2:6+i*3,1)=vP3D(:,i);

end

data{1}=pts1_n;
data{2}=pts2_n;

% options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt',...
%  'Display','iter','MaxIter',50);
% 
% [a,resnorm]=lsqnonlin(@loss_R_T,a0,[],[],options,data);



fprintf('global optimization:\n');
TolX=1e-6;
TolFun=1e-6;
MaxIter=40;
ConstantValue=[];
[a,resnorm]=Optimize_my_GN(@loss_R_T,@plus_function,a0,data,TolX,TolFun,MaxIter,ConstantValue);



phi=a(1:3);
% R1=eye(3);
% T1=[0;0;0];
R=angleAxisToRotationMatrix(phi);
T=a(4:6);


end



function [E,Jacbi]=loss_R_T(a,data)

pts1_n=data{1};
pts2_n=data{2};

n_P=size(pts1_n,1);

E=zeros(n_P*4+1,1);

Jacbi=zeros(n_P*4+1,size(a,1));

angleAxis2=a(1:3);
R1=eye(3);
T1=[0;0;0];
%R2=angleAxisToRotationMatrix(angleAxis2);
Q2=angleAxis2Quaternion(angleAxis2);
R2=quatern2rotMat(Q2);
T2=a(4:6);

for i=1:n_P

    P=a(6+i*3-2:6+i*3);

    V1=R1'*(P-T1);
    V2=R2'*(P-T2);

    E1=V1(1:2)/V1(3)-pts1_n(i,:)';

    E2=V2(1:2)/V2(3)-pts2_n(i,:)';

    E(i*4-3:i*4)=[E1;E2];

    d_E1_d_V1=[1/V1(3),0,-V1(1)/V1(3)^2;...
               0,1/V1(3),-V1(2)/V1(3)^2];

    d_E2_d_V2=[1/V2(3),0,-V2(1)/V2(3)^2;...
               0,1/V2(3),-V2(2)/V2(3)^2];

    d_V1_d_P=R1';

    d_V2_d_delta_theta=D_Rtrans_a_D_delta_theta(Q2,P-T2);

    d_V2_d_T2=-R2';

    d_V2_d_P=R2';

    Jacbi(i*4-3:i*4-2,6+i*3-2:6+i*3)=d_E1_d_V1*d_V1_d_P;


    Jacbi(i*4-1:i*4,6+i*3-2:6+i*3)=d_E2_d_V2*d_V2_d_P;


    Jacbi(i*4-1:i*4,1:3)=d_E2_d_V2*d_V2_d_delta_theta;

    Jacbi(i*4-1:i*4,4:6)=d_E2_d_V2*d_V2_d_T2;

end

norm_T2=norm(T2);

E(n_P*4+1)=norm_T2-1;

Jacbi(n_P*4+1,4:6)=1/norm_T2*[T2(1),T2(2),T2(3)];

end



function a=plus_function(a,delta_a,data)

%size_point=size(data{1},1);
if norm(delta_a(1:3,1))~=0

    Q1=quaternProd(angleAxis2Quaternion(a(1:3,1)), utility_deltaQ_VINS_Mono(delta_a(1:3,1)) );
    angleAxis1=quaternionToAngleAxis(Q1);
else

    angleAxis1=a(1:3,1);
end

a(1:3,1)=angleAxis1;

a(4:end,1)=a(4:end,1)+delta_a(4:end,1);
    
end








