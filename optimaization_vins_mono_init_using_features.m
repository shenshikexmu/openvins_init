function [x_I_k_opti,G_p_f_opti]=optimaization_vins_mono_init_using_features(x_I_k,features,camK,camD,map_camera_times)


size_frame=size(x_I_k,2);



for i=1:size_frame

    q_I_k=x_I_k(1:4,i);

    angleAxis_I_k=quaternionToAngleAxis(q_I_k);       % q_I_k_=angleAxis2Quaternion(angleAxis_I_k)

    a0(i*6-5:i*6)=[angleAxis_I_k;x_I_k(4:6,i)];

end

validateFeature=true(size(features,1),1);

all_ids = keys(features);

for j = 1:length(all_ids)

    id = all_ids{j}; % 获取当前键
  
    feat = features(id);
    
    featid=feat.featid;

    A=[];
    b=[];

    for frame=1:size_frame

        timestamp=map_camera_times(frame,1)-map_camera_times(frame,3);

        for cam_id = 1:length(feat.timestamps)

            for s=1:size(feat.timestamps{cam_id},1)

                if feat.timestamps{cam_id}(s,1)==timestamp

                    R=quatern2rotMat(x_I_k(1:4,frame));

                    T=x_I_k(5:7,frame);

                    uv_n=feat.uvs_norm{cam_id}(s,:);

                    V=[uv_n(1);uv_n(2);1];

                    A=[A;Skew_symmetric(R*V)];

                    b=[b;Skew_symmetric(R*V)*T];

                end

            end

        end

    end

    X=A\b;

    G_p_f(:,j)=X;

    for frame=1:size_frame

        timestamp=map_camera_times(frame,1)-map_camera_times(frame,3);

        for cam_id = 1:length(feat.timestamps)
        
            for s=1:size(feat.timestamps{cam_id},1)
            
                if feat.timestamps{cam_id}(s,1)==timestamp

                    R=quatern2rotMat(x_I_k(1:4,frame));
            
                    T=x_I_k(5:7,frame);
            
                    V=R'*(X-T);
                    if V(3)<0

                        validateFeature(j)=false;
                    end

                end

            end
        end
  
    end

end



size_ids=sum(validateFeature);

a0=zeros(size_frame*6+size_ids*3,1);



for j=1:size_ids

    a0(size_frame*6+j*3-2:size_frame*6+j*3)=G_p_f(:,j);

end






prior_Info=eye(10);

prior_Info(1:4,1:4)=prior_Info(1:4,1:4)*1/(1e-5)^2;           % 4dof unobservable yaw and position
prior_Info(5:7,5:7)=prior_Info(5:7,5:7)*1/(0.05)^2;           % bias_g prior
prior_Info(8:10,8:10)=prior_Info(8:10,8:10)*1/(0.10)^2;       % bias_a prior

prior_grad=zeros(10,1);

x_lin=x_I_k(:,1);


data{1}=imuPropagate;
data{2}=features;
data{3}=map_camera_times;
data{4}=camK;
data{5}=camD;
data{6}=camR;
data{7}=camT;
data{8}=size_frame;
data{9}=size_ids;
data{10}=[0;0;1]*gravity_mag;
data{11}=num_measurements;
data{12}=prior_Info;
data{13}=prior_grad;
data{14}=x_lin;
% data{10}=x_I_k;
% data{11}=G_p_f;


% options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',20);
% 
% [a,resnorm]=lsqnonlin(@loss_function_features,a0,[],[],options,data);
% 


fprintf('global optimization:\n');
TolX=1e-6;
TolFun=1e-6;
MaxIter=100;
ConstantValue=[4,5,6];
[a,resnorm]=Optimize_my_LM2(@loss_function_features,@plus_function,a0,data,TolX,TolFun,MaxIter,ConstantValue);



for j=1:size_ids

    G_p_f_opti(:,j)=a(size_frame*15+j*3-2:size_frame*15+j*3);

end


for n=1:size_frame

   a1=a(n*15-14:n*15,1);

   q_I_k_=angleAxis2Quaternion(a1(1:3));

   x_I_k_opti(:,n)=[q_I_k_;a1(4:15)];


end


end




function [E,Jacbi]=loss_function_features(a,data)
    
    
 

imuPropagate=data{1};
features=data{2};
map_camera_times=data{3};
camK=data{4};
camD=data{5};
camR=data{6};
camT=data{7};
size_frame=data{8};
size_ids=data{9};
gravity=data{10};
num_measurements=data{11};
prior_Info=data{12};
prior_grad=data{13};
x_lin=data{14};


E=zeros(num_measurements*2+(size_frame-1)*15+size(prior_grad,1),1);
Jacbi=zeros(num_measurements*2+(size_frame-1)*15+size(prior_grad,1),size_frame*15+size_ids*3);

%  Factor_ImuCPIv1
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
    

% Factor_ImageReprojCalib
precision_pic=1;
num_measurements=0;

all_ids = keys(features);

for i = 1:length(all_ids)
    id = all_ids{i}; % 获取当前键
  
    feat = features(id);
    
    featid=feat.featid;
    
    G_p_f_k=a(size_frame*15+(i-1)*3+1:size_frame*15+i*3);
    
    for cam_id = 1:length(feat.uvs_norm)
        uvs_norm=feat.uvs_norm{cam_id};
        uvs=feat.uvs{cam_id};
        timestamps=feat.timestamps{cam_id};
        
        for j=1:size(timestamps,1)
            
            num_measurements=num_measurements+1;
     
            camtime=timestamps(j,1);  
            n=find_cam_n_from_map_camera_times(camtime,map_camera_times);
            uv_norm=uvs_norm(j,:);
            uv=uvs(j,:);
            
            angleAxis1=a((n-1)*15+1:(n-1)*15+3,1);
            P1=a((n-1)*15+4:(n-1)*15+6,1);

            if n==1
                P1=[0;0;0];
            end
            
            Q1=angleAxis2Quaternion(angleAxis1);
            
            [Erepro,H_dz_dP1,H_dz_dangleAxis,H_dz_dG_p_f_k]=evaluate_Reprojection(P1,Q1,G_p_f_k,camK,camD,camR,camT,uv);
            
            E((size_frame-1)*15+num_measurements*2-1:(size_frame-1)*15+num_measurements*2,1)=Erepro*precision_pic;
            
            Jacbi((size_frame-1)*15+num_measurements*2-1:(size_frame-1)*15+num_measurements*2,(n-1)*15+1:(n-1)*15+3)=H_dz_dangleAxis*precision_pic;
            
            Jacbi((size_frame-1)*15+num_measurements*2-1:(size_frame-1)*15+num_measurements*2,(n-1)*15+4:(n-1)*15+6)=H_dz_dP1*precision_pic;
            
            Jacbi((size_frame-1)*15+num_measurements*2-1:(size_frame-1)*15+num_measurements*2,size_frame*15+(i-1)*3+1:size_frame*15+i*3)=H_dz_dG_p_f_k*precision_pic;
                    
        end
        
    end
    
    
end

% Factor_GenericPrior
% Comes from the form: cost = A * (x - x_lin) + b
llt0fI = chol(prior_Info, 'lower');
sqrtI=llt0fI';
b=sqrtI*prior_grad;

a1=a(1:15,1);
P1=a1(4:6,1);
angleAxis1=a1(1:3,1);
Q1=angleAxis2Quaternion(angleAxis1)';
%V1=a1(7:9,1);
Bg1=a1(10:12,1);
Ba1=a1(13:15,1);

Q_lin=x_lin(1:4,1);
ez=[0;0;1];

theta_err=quaternionToAngleAxis(quaternProd(invQuaternion(Q_lin),Q1));

res_Q1=sqrtI(1,1)*ez'*theta_err+b(1);

dq1_dtheta1=LeftMultiply( quaternProd(invQuaternion(Q_lin) , Q1)  );

J_Q1=sqrtI(1,1)*ez'*dq1_dtheta1(2:4,2:4);

res_P1=sqrtI(2:4,2:4)*(P1-x_lin(5:7,1))+b(2:4);
J_P1=sqrtI(2:4,2:4)*eye(3);

res_Bg=sqrtI(5:7,5:7)*(Bg1-x_lin(11:13,1))+b(5:7);
J_Bg=sqrtI(5:7,5:7)*eye(3);

res_Ba=sqrtI(8:10,8:10)*(Ba1-x_lin(14:16,1))+b(8:10);
J_Ba=sqrtI(8:10,8:10 )*eye(3);

E(num_measurements*2+(size_frame-1)*15+1:num_measurements*2+(size_frame-1)*15+size(prior_grad,1),1)=[res_Q1;res_P1;res_Bg;res_Ba] ;

Jacbi(num_measurements*2+(size_frame-1)*15+1:num_measurements*2+(size_frame-1)*15+1,1:3)=J_Q1;
Jacbi(num_measurements*2+(size_frame-1)*15+2:num_measurements*2+(size_frame-1)*15+4,4:6)=J_P1;
Jacbi(num_measurements*2+(size_frame-1)*15+5:num_measurements*2+(size_frame-1)*15+7,10:12)=J_Bg;
Jacbi(num_measurements*2+(size_frame-1)*15+8:num_measurements*2+(size_frame-1)*15+10,13:15)=J_Ba;




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










