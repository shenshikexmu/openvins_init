function [R,T]=optimaization_R_T(R,T,pts1_n,pts2_n,pts1,pts2,camK,camD)


q=rotMat2qRichard(R);

angleAxis=quaternionToAngleAxis(q);

%R_=angleAxisToRotationMatrix(angleAxis);


a0=[angleAxis;T];

R1=eye(3);
T1=[0;0;0];
R2=R;
T2=T;

validdata=true(size(pts1_n,1),1);

for i=1:size(pts1_n,1)


    U1=[pts1_n(i,1);pts1_n(i,2);1];
    
    U2=[pts2_n(i,1);pts2_n(i,2);1];
    
%     A=[Skew_symmetric(R1*U1);Skew_symmetric(R2*U2)];
% 
%     b=[Skew_symmetric(R1*U1)*T1;Skew_symmetric(R2*U2)*T2];
% 
%     P=A\b;

    A=[R1',-U1,U1*0;...
       R2',U2*0,-U2];

    b=[R1'*T1;R2'*T2];

    X=A\b;

    if X(4)>0 && X(5)>0

        P=X(1:3);
        a0=[a0;P];
    else

        validdata(i)=false;

    end

end


pts1_n=pts1_n(validdata,:);
pts2_n=pts2_n(validdata,:);
pts1=pts1(validdata,:);
pts2=pts2(validdata,:);


data{1}=pts1_n;
data{2}=pts2_n;
data{3}=pts1;
data{4}=pts2;
data{5}=camK;
data{6}=camD;


options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',20);

[a,resnorm]=lsqnonlin(@loss_R_T,a0,[],[],options,data);


angleAxis=a(1:3);

R=angleAxisToRotationMatrix(angleAxis);

T=a(4:6);

if 0

    figure

    colors = jet(size(R,2));
    
    RR{1}=eye(3);
    TT{1}=zeros(3,1);
    RR{2}=R;
    TT{2}=T;
    
    scale=0.5;
    
    img_x_len=0.2*scale;
    img_y_len=0.127*scale;
    
    for i=1:size(RR,2)
    
        G_camR_k=RR{i};  
        G_camT_k=TT{i};
        
        p1=G_camT_k-G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
        p2=G_camT_k-G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
        p3=G_camT_k+G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
        p4=G_camT_k+G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
        
        lines=[p1,p2,p3,p4,p1];
        
        plot3(lines(1,:),lines(2,:),lines(3,:),'-', 'Color', colors(i,:));
        hold on;
        if i==1
    
            text(G_camT_k(1),G_camT_k(2),G_camT_k(3),'0');
            hold on
    
            plot3([0,1]*0.1*scale,[0,0],[0,0],'r-',...
                  [0,0],[0,1]*0.1*scale,[0,0],'g-',...
                  [0,0],[0,0],[0,1]*0.1*scale,'b-' );
    
        end
         
    end
    
    axis equal;
    
    Len=16;
    
    for i = 1:size(pts1_n,1)
    
    
        V1=RR{1}*[pts1_n(i,1);pts1_n(i,2);1];
        V2=RR{2}*[pts2_n(i,1);pts2_n(i,2);1];
    
        P=a(6+i*3-2:6+i*3);
                       
        plot3([TT{1}(1),TT{1}(1)+V1(1)*Len],[TT{1}(2),TT{1}(2)+V1(2)*Len],[TT{1}(3),TT{1}(3)+V1(3)*Len],'r-',...
              [TT{2}(1),TT{2}(1)+V2(1)*Len],[TT{2}(2),TT{2}(2)+V2(2)*Len],[TT{2}(3),TT{2}(3)+V2(3)*Len],'g-',...
              P(1),P(2),P(3),'b.');
    
    end
    
    hold off;

end


end





function E=loss_R_T(a,data)

pts1_n=data{1};
pts2_n=data{2};

pts1=data{3};
pts2=data{4};
camK=data{5};
camD=data{6};

angleAxis=a(1:3);

R=angleAxisToRotationMatrix(angleAxis);

T=a(4:6);

R1=eye(3);
T1=[0;0;0];
R2=R;
T2=T;

E=zeros(size(pts1_n,1)*4+1,1);

for i=1:size(pts1_n,1)

    P=a(6+i*3-2:6+i*3);

    V1=R1'*(P-T1);
    V2=R2'*(P-T2);

%    E(i*6-5:i*6)=[Skew_symmetric(V1)*[pts1_n(i,1:2),1]';Skew_symmetric(V2)*[pts2_n(i,1:2),1]'];

%    E(i*4-3:i*4)=[V1(1:2)/V1(3)-pts1_n(i,1:2)';V2(1:2)/V2(3)-pts2_n(i,1:2)'];

    [uv1_dist,~]= distort_cv(V1(1:2)/V1(3), camK,camD);
    [uv2_dist,~]= distort_cv(V2(1:2)/V2(3), camK,camD);

    E1=uv1_dist-pts1(i,1:2);
    E2=uv2_dist-pts2(i,1:2);

    if V1(3)<0
        E1=[1000,1000];
    end

    if V2(3)<0
        E2=[1000,1000];
    end

    E(i*4-3:i*4)=[E1,E2]';

end

E(size(pts1_n,1)*4+1)=(norm(T)-1)*1;


end







