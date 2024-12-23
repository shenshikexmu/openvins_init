function [R,T]=Initial_R_T(P1,P2)


PP1=P1(:,1:2);
PP2=P2(:,1:2);

n=size(PP1,1);


A=[];

for i=1:n
    
    temp_A=[PP1(i,1)*PP2(i,1), PP1(i,1)*PP2(i,2), PP1(i,1), PP1(i,2)*PP2(i,1),...
            PP1(i,2)*PP2(i,2), PP1(i,2), PP2(i,1), PP2(i,2), 1];
    
    A=[A;temp_A];
    
end

%[UU,SS,VV] = svd(A);

%[~,~,VV2] = svd(A,'econ');

[UU2,SS2,VV2] = svd(A);

e=VV2(:,end);

E=[ e(1),e(2),e(3);...
    e(4),e(5),e(6);...
    e(7),e(8),e(9)];

[U,S1,V]=svd(E);
%[U2,S2,V2]=svd(-E);


if det(U)<-0.95
    U=-U;
end
if det(V)<-0.95
    V=-V;
end

Rz90=[0,-1,0;...
      1, 0,0;...
      0, 0,1];




T1=U(:,3);
R1=U*Rz90*V';
T2=-U(:,3);
R2=U*Rz90*V';
T3=U(:,3);
R3=U*Rz90'*V';
T4=-U(:,3);
R4=U*Rz90'*V';



[R,T]=Judge_right_R_T2(PP1(1,1:2),PP2(1,1:2),R1,R2,R3,R4,T1,T2,T3,T4);





end


function [R,T]=Judge_right_R_T2(P1,P2,R1,R2,R3,R4,T1,T2,T3,T4)

R=[];
T=[];


[s11,s12]=Location_scale(P1,P2,R1,T1);
[s21,s22]=Location_scale(P1,P2,R2,T2);
[s31,s32]=Location_scale(P1,P2,R3,T3);
[s41,s42]=Location_scale(P1,P2,R4,T4);
%SS=[s11,s12;s21,s22;s31,s32;s41,s42];
if s11>0 && s12>0
     R=R1;
     T=T1;
end
if s21>0 && s22>0
     R=R2;
     T=T2;
end
if s31>0 && s32>0
     R=R3;
     T=T3;
end
if s41>0 && s42>0
     R=R4;
     T=T4;
end



end



function [s1,s2]=Location_scale(P1,P2,R,T)

V1=[P1,1]';
V2=R*[P2,1]';

t=[V1,-V2]\T;

s1=t(1);
s2=t(2);

end



function Skew=Skew_symmetric(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end





