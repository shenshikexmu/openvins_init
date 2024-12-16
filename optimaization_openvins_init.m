function [x_I_k_opti,G_p_f_opti]=optimaization_openvins_init(x_I_k,G_p_f,imuPropagate,pts_opti,pts_n_opti,camK,camD,camR,camT,gravity_mag)


size_frame=size(x_I_k,2);

size_ids=size(G_p_f,2);

a0=zeros(size_frame*15+size_ids*3,1);

for i=1:size_frame

    q_I_k=x_I_k(1:4,i);

    angleAxis_I_k=quaternionToAngleAxis(q_I_k);       % q_I_k_=angleAxis2Quaternion(angleAxis_I_k)

    a0(i*15-14:i*15)=[angleAxis_I_k;x_I_k(5:16,i)];

end


for j=1:size_ids

    a0(size_frame*15+j*3-2:size_frame*15+j*3)=G_p_f(:,j);

end


data{1}=imuPropagate;
data{2}=pts_opti;
data{3}=pts_n_opti;
data{4}=camK;
data{5}=camD;
data{6}=camR;
data{7}=camT;
data{8}=size_frame;
data{9}=size_ids;
data{10}=[0;0;1]*gravity_mag;
% data{10}=x_I_k;
% data{11}=G_p_f;


% options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',20);
% 
% [a,resnorm]=lsqnonlin(@loss_function,a0,[],[],options,data);
% 


fprintf('global optimization:\n');
TolX=1e-6;
TolFun=1e-6;
MaxIter=100;
ConstantValue=[4,5,6];
[a,resnorm]=Optimize_my_LM2(@loss_function,@plus_function,a0,data,TolX,TolFun,MaxIter,ConstantValue);



for j=1:size_ids

    G_p_f_opti(:,j)=a(size_frame*15+j*3-2:size_frame*15+j*3);

end


for n=1:size_frame

   a1=a(n*15-14:n*15,1);


   q_I_k_=angleAxis2Quaternion(a1(1:3));

   x_I_k_opti(:,n)=[q_I_k_;a1(4:15)];


end


end




function [E,Jacbi]=loss_function(a,data)

imuPropagate=data{1};
pts_opti=data{2};
pts_n_opti=data{3};
camK=data{4};
camD=data{5};
camR=data{6};
camT=data{7};
size_frame=data{8};
size_ids=data{9};
gravity=data{10};


E=zeros(size_frame*size_ids*2+(size_frame-1)*15,1);
Jacbi=zeros(size_frame*size_ids*2+(size_frame-1)*15,size_frame*15+size_ids*3);


for i=2:size_frame
    
    %a0(i*15+1:i*15+15,1)=[angleAxistemp;Ptemp;Vtemp;Batemp;Bgtemp];
    
    a1=a((i-2)*15+1:(i-1)*15,1);%a(i*15-29:i*15-15,1);
    a2=a((i-1)*15+1:i*15,1);%a(i*15-14:i*15,1);
    
    P1=a1(4:6,1);
    angleAxis1=a1(1:3,1);
    Q1=angleAxis2Quaternion(angleAxis1)';
    V1=a1(7:9,1);
    Bg1=a1(10:12,1);
    Ba1=a1(13:15,1);
    
    P2=a2(4:6,1);
    angleAxis2=a2(1:3,1);
    Q2=angleAxis2Quaternion(angleAxis2)';
    V2=a2(7:9,1);
    Bg2=a2(10:12,1);
    Ba2=a2(13:15,1);
    
    [residuals,J1,J2]=evaluate_VINS_Mono_gravity(P1,Q1,V1,Ba1,Bg1,P2,Q2,V2,Ba2,Bg2,imuPropagate{i},gravity); 
    
    E((i-2)*15+1:(i-1)*15,1)=residuals;
    
    Jacbi((i-2)*15+1:(i-1)*15,(i-2)*15+1:(i-1)*15)=[J1(:,4:6),J1(:,1:3),J1(:,7:9),J1(:,13:15),J1(:,10:12)];
    
    Jacbi((i-2)*15+1:(i-1)*15,(i-1)*15+1:i*15)=[J2(:,4:6),J2(:,1:3),J2(:,7:9),J2(:,13:15),J2(:,10:12)];
 
end
    

precision_pic=1;


for j=1:size_ids

    G_p_f_k=a(size_frame*15+(j-1)*3+1:size_frame*15+j*3);

    for n=1:size_frame

        pts_temp=pts_opti{n}(j,1:2);

        pts_n_temp=pts_n_opti{n}(j,1:2);

%        a1=a((n-1)*15+1:(n-1)*15+15,1);
%        P1=a1(4:6,1);       
%        angleAxis1=a1(1:3,1);
        
        angleAxis1=a((n-1)*15+1:(n-1)*15+3,1);
        P1=a((n-1)*15+4:(n-1)*15+6,1);

        if n==1
            P1=[0;0;0];
        end
        
        Q1=angleAxis2Quaternion(angleAxis1);
        %R1=angleAxisToRotationMatrix(angleAxis1);     %  Q1=angleAxis2Quaternion(angleAxis1)';   R1_=quatern2rotMat(Q1)
        
        [Erepro,H_dz_dP1,H_dz_dangleAxis,H_dz_dG_p_f_k]=evaluate_Reprojection(P1,Q1,G_p_f_k,camK,camD,camR,camT,pts_temp);

        %Erepro=evaluate_Reprojection_norm(P1,R1,camK,camD,camR,camT,pts_n_temp,G_p_f_k)*precision_pic*camK(1,1);
  
        E((size_frame-1)*15+(j-1)*size_frame*2+(n-1)*2+1:(size_frame-1)*15+(j-1)*size_frame*2+n*2,1)=Erepro*precision_pic;
        
        Jacbi((size_frame-1)*15+(j-1)*size_frame*2+(n-1)*2+1:(size_frame-1)*15+(j-1)*size_frame*2+n*2,(n-1)*15+1:(n-1)*15+3)=H_dz_dangleAxis*precision_pic;
         
        Jacbi((size_frame-1)*15+(j-1)*size_frame*2+(n-1)*2+1:(size_frame-1)*15+(j-1)*size_frame*2+n*2,(n-1)*15+4:(n-1)*15+6)=H_dz_dP1*precision_pic;
        
        Jacbi((size_frame-1)*15+(j-1)*size_frame*2+(n-1)*2+1:(size_frame-1)*15+(j-1)*size_frame*2+n*2,size_frame*15+(j-1)*3+1:size_frame*15+j*3)=H_dz_dG_p_f_k*precision_pic;
        


    end


end




end





function a=plus_function(a,delta_a,data)


size_frame=data{8};
size_ids=data{9};


for i=1:size_frame
    
    %a0(i*15+1:i*15+15,1)=[angleAxistemp;Ptemp;Vtemp;Batemp;Bgtemp];
    
    a1=a((i-1)*15+1:i*15,1);
    delta_a1=delta_a((i-1)*15+1:i*15,1);

    if norm(delta_a1(1:3,1))~=0

        Q1=quaternProd(angleAxis2Quaternion(a1(1:3,1)), utility_deltaQ_VINS_Mono(delta_a1(1:3,1)) );
        angleAxis1=quaternionToAngleAxis(Q1);
    else

        angleAxis1=a1(1:3,1);
    end
    
    P1=a1(4:6,1)+delta_a1(4:6,1);

    V1=a1(7:9,1)+delta_a1(7:9,1);
    Bg1=a1(10:12,1)+delta_a1(10:12,1);
    Ba1=a1(13:15,1)+delta_a1(13:15,1);

    a((i-1)*15+1:i*15,1)=[angleAxis1;P1;V1;Bg1;Ba1];

   
end
    

for j=1:size_ids

    %G_p_f_k=a(size_frame*15+(j-1)*3+1:size_frame*15+j*3);

    a(size_frame*15+(j-1)*3+1:size_frame*15+j*3)=a(size_frame*15+(j-1)*3+1:size_frame*15+j*3)+delta_a(size_frame*15+(j-1)*3+1:size_frame*15+j*3);

end


end











function [E,H_dz_dG_I_p,H_dz_dG_I_delta_theta,H_dz_dG_p_f_k]=evaluate_Reprojection(G_I_p,G_I_q,G_p_f_k,camK,camD,camR,camT,pts_temp)

G_I_R=quatern2rotMat(G_I_q);

G_c_R=G_I_R*camR;

G_c_T=G_I_p+G_I_R*camT;

C_p_f=G_c_R'*(G_p_f_k-G_c_T);

[uv_dist,H_dz_dzn]= distort_cv(C_p_f(1:2)/C_p_f(3), camK,camD);

E=(uv_dist-pts_temp)';

d_uv_norm_d_C_p_f=[1/C_p_f(3),0,-C_p_f(1)/C_p_f(3)^2;...
                   0,1/C_p_f(3),-C_p_f(2)/C_p_f(3)^2];
                   
d_C_p_f_d_G_p_f_k=G_c_R';

d_C_p_f_d_G_I_p=-G_c_R';

d_C_p_f_d_G_I_delta_theta=camR'*D_Rtrans_a_D_delta_theta(G_I_q,G_p_f_k-G_I_p); % 问题在于 G_c_T包含 G_I_R，所以在计算导数时需要把G_c_T中的G_I_R考虑在内。 

H_dz_dG_p_f_k=H_dz_dzn*d_uv_norm_d_C_p_f*d_C_p_f_d_G_p_f_k;

H_dz_dG_I_p=H_dz_dzn*d_uv_norm_d_C_p_f*d_C_p_f_d_G_I_p;


H_dz_dG_I_delta_theta=H_dz_dzn*d_uv_norm_d_C_p_f*d_C_p_f_d_G_I_delta_theta;


end


function E=evaluate_Reprojection_norm(G_I_p,G_I_R,camK,camD,camR,camT,pts_n_temp,G_p_f_k)



G_c_R=G_I_R*camR;

G_c_T=G_I_p+G_I_R*camT;


C_p_f=G_c_R'*(G_p_f_k-G_c_T);

C_p_f=C_p_f/C_p_f(3);


Gamma=[1,0,-pts_n_temp(1);...
       0,1,-pts_n_temp(2)];


E=Gamma*C_p_f;


end










