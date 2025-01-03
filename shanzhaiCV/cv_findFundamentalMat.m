function [mask,R,T]=cv_findFundamentalMat(points1, points2, method, ransacReprojThreshold ,confidence )
    
n=size(points1,1);
mask=zeros(n,1);


if  strcmp(method , 'cv_FM_RANSAC')


    n_linear=min(8,max(8,floor(n/2)));            %需要解方程的行数

    k=ceil(log(1-confidence)/log(1-0.8^n_linear));   %需要循环的次数

    k=200;

   
%    min_SS=inf;

    max_SS=[];
    max_count=0;

    min_SS=inf;

    for i=1:k
        
        idx = randperm(size(points1,1),n_linear); 
        P1_temp= points1(idx,1:2);
        P2_temp= points2(idx,1:2);
        [R,T]=Initial_R_T(P1_temp,P2_temp);
        SS=zeros(size(points1,1),1);
        flag=zeros(size(points1,1),1);
        if ~isempty(R)
            T=T/norm(T);
            EEE=Skew_symmetric(T)*R;
            for s=1:size(points1,1)
                SS(s,1)=computeError(points1(s,1:2),points2(s,1:2),EEE);%abs([points1(s,1:2),1]*EEE*[points2(s,1:2),1]');             
                [s1,s2]=Location_scale(points1(s,1:2),points2(s,1:2),R,T);
                if s1>0 &&s2>0
                    flag(s,1)=1;
                end

            end
    
            count = sum(flag);
            

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

    EEE=Skew_symmetric(T)*R;


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
%     %m=ceil(n*0.8/12);
% 
%     P1_temp= points1(sort_idx(1:floor(n/2)),1:2);
%     P2_temp= points2(sort_idx(1:floor(n/2)),1:2);
%     [R,T]=Initial_R_T(P1_temp,P2_temp);
% 
%     EEE=Skew_symmetric(T)*R;
% 
% 
% 
%     for s=1:size(points1,1)
% 
%          SS(s,1)=computeError(points1(s,1:2),points2(s,1:2),EEE);
% 
% 
%         if (computeError(points1(s,1:2),points2(s,1:2),EEE)<ransacReprojThreshold)
% 
%             mask(s,1)=1;
% 
%         else
% 
%             SS(s,1)=0;
% 
%         end
% 
%     end




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
    
    
    %E=Skew_symmetric(T)*R;
    

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


function [s1,s2]=Location_scale(P1,P2,R,T)

V1=[P1,1]';
V2=R*[P2,1]';

t=[V1,-V2]\T;

s1=t(1);
s2=t(2);

end

