function [R,T]=Initial_R_T_CameraHouse_Tsing_RANSAC_fisheye3(PP1,FishIntrinsic1,PP2,FishIntrinsic2,L_light,inlier_rate)

P1=PP_distortion_fisheye(PP1,FishIntrinsic1);

P2=PP_distortion_fisheye(PP2,FishIntrinsic2);

n_linear=14;            %需要解方程的行数

k=ceil(log(1-0.99999)/log(1-inlier_rate^n_linear));   %需要循环的次数

n_r=0;
for i=1:size(P1,1)
    if P1(i,2)~=0 && P2(i,2)~=0
        n_r=n_r+1;
        PPr1(n_r,1:2)=P1(i,1:2);
        PPr2(n_r,1:2)=P2(i,1:2);
    end
end

n_b=0;
for i=1:size(P1,1)
    if P1(i,4)~=0 && P2(i,4)~=0
        n_b=n_b+1;
        PPb1(n_b,1:2)=P1(i,3:4);
        PPb2(n_b,1:2)=P2(i,3:4);
    end
end



%max_SS=0;
min_SS=inf;

for i=1:k
    
    idx = randperm(n_r,n_linear); 
    P1_temp= PPr1(idx,1:2);
    P2_temp= PPr2(idx,1:2);
    [R,T]=Initial_R_T_CameraHouse2(P1_temp,P2_temp);
    if ~isempty(R)
        T=T/norm(T);
        EEE=V_2_Skew(T)*R;
        for s=1:size(PPb1,1)
            SS(s,1)=abs([PPb1(s,1:2),1]*EEE*[PPb2(s,1:2),1]');

        end

%         total=sum(SS<0.03);
%         if total>max_SS
%             maz_SS=total;
%             max_R=R;
%             max_T=T;
%         end
        
        sum_SS=sum(SS);
        if sum_SS<min_SS
            min_SS=sum_SS;
            min_R=R;
            min_T=T;
        end
    end
end



R=min_R;
T=min_T;



EEE=V_2_Skew(T)*R;

nn=0;

for i=1:size(P1,1)
    
     if P1(i,2)~=0 && P2(i,2)~=0 && P1(i,4)~=0 && P2(i,4)~=0
         
         nn=nn+1;
         Xr1=Location_X(P1(i,1:2),P2(i,1:2),R,T);   
         Xb1=Location_X(P1(i,3:4),P2(i,3:4),R,T);
         L(nn,1:2)=[norm(Xr1-Xb1),i];
         
     end
    
end

L_median=median(L(:,1));



scale=L_light/L_median;

T=T*scale;


end

function PP_distort=PP_distortion_fisheye(PP,FishIntrinsic)

for i=1:size(PP,1)
  
    if PP(i,2)~=0

        PP_distort(i,1:2) = ...
            My_undistortFisheyePoints(PP(i,1:2),FishIntrinsic);

    else
            PP_distort(i,1:2)=[0,0];
      
    end
    
   if PP(i,4)~=0

        PP_distort(i,3:4) = ...
            My_undistortFisheyePoints(PP(i,3:4),FishIntrinsic);

    else
            PP_distort(i,3:4)=[0,0];
      
    end


end
end


function Skew=V_2_Skew(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end




function [X,s1,s2]=Location_X(P1,P2,R,T)

V1=[P1,1]';
V2=R*[P2,1]';

A=[V1(2), -V1(1),   0;...
   V1(3),   0,    -V1(1);...
    0,    V1(3),  -V1(2);...
   V2(2), -V2(1),   0;...
   V2(3),   0,    -V2(1);...
    0,    V2(3),  -V2(2)];

B=[0;...
   0;...
   0;...
   V2(2)*T(1)-V2(1)*T(2);...
   V2(3)*T(1)-V2(1)*T(3);...
   V2(3)*T(2)-V2(2)*T(3)];

X=A\B;


s1=X(3);



s2=(X-T)'*R(:,3);

end


function [R,T]=Initial_R_T_CameraHouse2(P1,P2)


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
S=[1,0,0;...
   0,1,0;...
   0,0,0 ];




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
SS=[s11,s12;s21,s22;s31,s32;s41,s42];
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





