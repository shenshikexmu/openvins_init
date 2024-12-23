function [R,T]=optimaization_R_T(R,T,PP1,PP2)


q=rotMat2qRichard(R);

angleAxis=quaternionToAngleAxis(q);

R_=angleAxisToRotationMatrix(angleAxis);


a0=[angleAxis;T];


data{1}=PP1;
data{2}=PP2;


options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',20);

[a,resnorm]=lsqnonlin(@loss_R_T,a0,[],[],options,data);

angleAxis=a(1:3);

R=angleAxisToRotationMatrix(angleAxis);

T=a(4:6);



end


function E=loss_R_T(a,data)

PP1=data{1};
PP2=data{2};

angleAxis=a(1:3);

R=angleAxisToRotationMatrix(angleAxis);

T=a(4:6);

E=zeros(size(PP1,1)+1,1);

for i=1:size(PP1,1)

    V1=[PP1(i,1);PP1(i,2);1];

    V2=[PP2(i,1);PP2(i,2);1];

    E(i)=V1'*Skew_symmetric(T)*R*V2;


end

E(size(PP1,1)+1)=(norm(T)-1)*1000;





end