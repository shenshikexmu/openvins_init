function [mask,R,T]=cv_findFundamentalMat(points1, points2, method, ransacReprojThreshold ,confidence )
    
n=size(points1,1);
mask=zeros(n,1);


if  strcmp(method , 'cv_FM_RANSAC')


    n_linear=min(8,max(8,floor(n/2)));            %需要解方程的行数

    k=ceil(log(1-confidence)/log(1-0.8^n_linear))   %需要循环的次数


%    min_SS=inf;

    max_SS=[];
    max_count=0;

    min_SS=inf;

    for i=1:k
        
        idx = randperm(size(points1,1),n_linear); 
        P1_temp= points1(idx,1:2);
        P2_temp= points2(idx,1:2);
        [R,T]=Initial_R_T(P1_temp,P2_temp);
        if ~isempty(R)
            T=T/norm(T);
            EEE=V_2_Skew(T)*R;
            for s=1:size(points1,1)
                SS(s,1)=computeError(points1(s,1:2),points2(s,1:2),EEE);%abs([points1(s,1:2),1]*EEE*[points2(s,1:2),1]');
    
            end
    
    %         total=sum(SS<0.03);
    %         if total>max_SS
    %             maz_SS=total;
    %             max_R=R;
    %             max_T=T;
    %         end
            
%             sum_SS=sum(SS);
% %                 i
% %                 sum_SS
%             if sum_SS<min_SS
%                 min_SS=sum_SS;
%                 min_R=R;
%                 min_T=T;
%             end

            count = sum(SS < ransacReprojThreshold);

            if count>max_count

                max_count=count;

                min_SS=sum(SS);

                max_R=R;
                max_T=T;


            elseif count==max_count

                sum_SS=sum(SS);

                if sum_SS<min_SS

                    min_SS=sum_SS;

                    max_R=R;
                    max_T=T;

                end


            end

        end
    end



%     R=min_R;
%     T=min_T;

    R=max_R;
    T=max_T;

    EEE=V_2_Skew(T)*R;


    for s=1:size(points1,1)

         SS(s,1)=computeError(points1(s,1:2),points2(s,1:2),EEE);


%         if (abs([points1(s,1:2),1]*EEE*[points2(s,1:2),1]')<ransacReprojThreshold)
% 
%             mask(s,1)=1;
% 
%         else
% 
%             SS(s,1)=0;
% 
%         end

    end


    [sortSS,sort_idx]=sort(SS);

    %m=ceil(n*0.8/12);

    P1_temp= points1(sort_idx(1:floor(n/2)),1:2);
    P2_temp= points2(sort_idx(1:floor(n/2)),1:2);
    [R,T]=Initial_R_T(P1_temp,P2_temp);

    EEE=V_2_Skew(T)*R;



    for s=1:size(points1,1)

         SS(s,1)=computeError(points1(s,1:2),points2(s,1:2),EEE);


        if (computeError(points1(s,1:2),points2(s,1:2),EEE)<ransacReprojThreshold)

            mask(s,1)=1;

        else

            SS(s,1)=0;

        end

    end




%     [sortSS,sort_idx]=sort(SS);
% 
%     for s=n:-1:n-4
% 
%         mask(sort_idx(s,1),1)=0;
% 
%         SS(sort_idx(s,1),1)=0;
% 
% 
%     end
% 
% 


    
end


    
    
    
end



function  err=computeError(p1,p2,E)
    
    
    %E=V_2_Skew(T)*R;
    

    v1=[p1(1);p1(2);1];
    v2=[p2(1);p2(2);1];
    
    d1=v1'*E*v2;
    
    v1_=E*v2;
    
    v2_=E'*v1;
    
    d2=d1;
    
    s1=1/(v1_(1)^2+v1_(2)^2);
    
    s2=1/(v2_(1)^2+v2_(2)^2);
    
    
    err=max(d1*d1*s1,d2*d2*s2);
    
    err=sqrt(err);
   
end







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



function Skew=V_2_Skew(V)

Skew=[   0,  -V(3),  V(2);...
       V(3),    0,  -V(1);...
      -V(2),  V(1),    0];


end





