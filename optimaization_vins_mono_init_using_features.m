function [x_I_k_opti,G_p_f_opti]=optimaization_vins_mono_init_using_features(x_I_k,features,camK,camD,map_camera_times)


size_frame=size(x_I_k,2);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     for frame=1:size_frame
% 
%         timestamp=map_camera_times(frame,1)-map_camera_times(frame,3);
% 
%         for cam_id = 1:length(feat.timestamps)
%         
%             for s=1:size(feat.timestamps{cam_id},1)
%             
%                 if feat.timestamps{cam_id}(s,1)==timestamp
% 
%                     R=quatern2rotMat(x_I_k(1:4,frame));
%             
%                     T=x_I_k(5:7,frame);
%             
%                     V=R'*(X-T);
%                     if V(3)<0
% 
%                         validateFeature(j)=false;
%                     end
% 
%                 end
% 
%             end
%         end
%   
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

for j = 1:length(all_ids)
    id = all_ids{j}; % 获取当前键
    if validateFeature(j)==false
        remove(features, id);
    end
end

size_ids=sum(validateFeature);

G_p_f=G_p_f(:,validateFeature);

a0=zeros(size_frame*6+size_ids*3,1);

for i=1:size_frame

    q_I_k=x_I_k(1:4,i);

    angleAxis_I_k=quaternionToAngleAxis(q_I_k);       % q_I_k_=angleAxis2Quaternion(angleAxis_I_k)

    a0(i*6-5:i*6)=[angleAxis_I_k;x_I_k(5:7,i)];

end

for j=1:size_ids

    a0(size_frame*6+j*3-2:size_frame*6+j*3)=G_p_f(:,j);

end

num_measurements=0;
all_ids = keys(features);

for j = 1:length(all_ids)

    id = all_ids{j}; % 获取当前键
  
    feat = features(id);

    for cam_id = 1:length(feat.timestamps)

        num_measurements=num_measurements+size(feat.timestamps{cam_id},1);

    end

end


prior_Info=eye(6);
prior_Info(1:6,1:6)=prior_Info(1:6,1:6)*1/(1e-5)^2;           % 6dof unobservable yaw and position
prior_grad=zeros(6,1);

x_lin=x_I_k(:,1);
camR=eye(3);
camT=zeros(3,1);

data{1}=[];
data{2}=features;
data{3}=map_camera_times;
data{4}=camK;
data{5}=camD;
data{6}=camR;
data{7}=camT;
data{8}=size_frame;
data{9}=size_ids;
data{10}=[];
data{11}=num_measurements;
data{12}=prior_Info;
data{13}=prior_grad;
data{14}=x_lin;
% data{10}=x_I_k;
% data{11}=G_p_f;


% options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',500);
% 
% [a,resnorm]=lsqnonlin(@loss_function_features,a0,[],[],options,data);



fprintf('global optimization:\n');
TolX=1e-6;
TolFun=1e-6;
MaxIter=200;
ConstantValue=[1,2,3,4,5,6];
[a,resnorm]=Optimize_my_LM2(@loss_function_features,@plus_function,a0,data,TolX,TolFun,MaxIter,ConstantValue);



for j=1:size_ids

    G_p_f_opti(:,j)=a(size_frame*6+j*3-2:size_frame*6+j*3);

end


for n=1:size_frame

   a1=a(n*6-5:n*6,1);

   q_I_k_opti=angleAxis2Quaternion(a1(1:3));

   x_I_k_opti(:,n)=[q_I_k_opti;a1(4:6)];

end


end




function [E,Jacbi]=loss_function_features(a,data)




features=data{2};
map_camera_times=data{3};
camK=data{4};
camD=data{5};
camR=data{6};
camT=data{7};
size_frame=data{8};
size_ids=data{9};
num_measurements=data{11};
prior_Info=data{12};
prior_grad=data{13};
x_lin=data{14};


E=zeros(num_measurements*2+size(prior_grad,1)+1,1);
Jacbi=zeros(num_measurements*2+size(prior_grad,1)+1,size_frame*6+size_ids*3);



% Factor_ImageReprojCalib
precision_pic=1;
num_measurements=0;

all_ids = keys(features);

for i = 1:length(all_ids)
    id = all_ids{i}; % 获取当前键
  
    feat = features(id);
    
    G_p_f_k=a(size_frame*6+(i-1)*3+1:size_frame*6+i*3);
    
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
            
            angleAxis1=a((n-1)*6+1:(n-1)*6+3,1);
            P1=a((n-1)*6+4:(n-1)*6+6,1);

%             if n==1
%                 angleAxis1=[0;0;0];
%                 P1=[0;0;0];
%             end
            
            Q1=angleAxis2Quaternion(angleAxis1);
            
            [Erepro,H_dz_dP1,H_dz_dangleAxis,H_dz_dG_p_f_k]=evaluate_Reprojection(P1,Q1,G_p_f_k,camK,camD,camR,camT,uv);

            %Erepro=evaluate_Reprojection_norm(P1,Q1,camR,camT,uv_norm,G_p_f_k)*500;
            
            E(num_measurements*2-1:num_measurements*2,1)=Erepro*precision_pic;
            
            Jacbi(num_measurements*2-1:num_measurements*2,(n-1)*6+1:(n-1)*6+3)=H_dz_dangleAxis*precision_pic;
            
            Jacbi(num_measurements*2-1:num_measurements*2,(n-1)*6+4:(n-1)*6+6)=H_dz_dP1*precision_pic;
            
            Jacbi(num_measurements*2-1:num_measurements*2,size_frame*6+(i-1)*3+1:size_frame*6+i*3)=H_dz_dG_p_f_k*precision_pic;
                    
        end
        
    end
    
    
end

% Factor_GenericPrior
% Comes from the form: cost = A * (x - x_lin) + b
llt0fI = chol(prior_Info, 'lower');
sqrtI=llt0fI';
b=sqrtI*prior_grad;

a1=a(1:6,1);

res_Q_P=(a1-zeros(6,1))+b;

E(num_measurements*2+1:num_measurements*2+size(prior_grad,1),1)=sqrtI(1,1)*res_Q_P;
Jacbi(num_measurements*2+1:num_measurements*2+size(prior_grad,1),1:6)=sqrtI(1,1)*eye(size(prior_grad,1));


P1=a(4:6,1);
P2=a(10:12,1);

norm_P1_P2=norm(P1-P2);

res_P1_P2=(norm_P1_P2-1);

d_norm_P1_P2_d_P1=1/norm_P1_P2*(P1-P2)';

d_norm_P1_P2_d_P2=-d_norm_P1_P2_d_P1;

E(num_measurements*2+size(prior_grad,1)+1:num_measurements*2+size(prior_grad,1)+1,1)=sqrtI(1,1)*res_P1_P2;

Jacbi(num_measurements*2+size(prior_grad,1)+1:num_measurements*2+size(prior_grad,1)+1,4:6)=sqrtI(1,1)*d_norm_P1_P2_d_P1;

Jacbi(num_measurements*2+size(prior_grad,1)+1:num_measurements*2+size(prior_grad,1)+1,10:12)=sqrtI(1,1)*d_norm_P1_P2_d_P2;



end





function a=plus_function(a,delta_a,data)


size_frame=data{8};
size_ids=data{9};


for i=1:size_frame
    
    %a0(i*15+1:i*15+15,1)=[angleAxistemp;Ptemp;Vtemp;Batemp;Bgtemp];
    
    a1=a((i-1)*6+1:i*6,1);
    delta_a1=delta_a((i-1)*6+1:i*6,1);

    if norm(delta_a1(1:3,1))~=0

        Q1=quaternProd(angleAxis2Quaternion(a1(1:3,1)), utility_deltaQ_VINS_Mono(delta_a1(1:3,1)) );
        angleAxis1=quaternionToAngleAxis(Q1);
    else

        angleAxis1=a1(1:3,1);
    end
    
    P1=a1(4:6,1)+delta_a1(4:6,1);

    a((i-1)*6+1:i*6,1)=[angleAxis1;P1];
 
end
    

for j=1:size_ids

    %G_p_f_k=a(size_frame*15+(j-1)*3+1:size_frame*15+j*3);

    a(size_frame*6+(j-1)*3+1:size_frame*6+j*3)=a(size_frame*6+(j-1)*3+1:size_frame*6+j*3)+delta_a(size_frame*6+(j-1)*3+1:size_frame*6+j*3);

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


function E=evaluate_Reprojection_norm(G_I_p,G_I_q,camR,camT,pts_n_temp,G_p_f_k)

G_I_R=quatern2rotMat(G_I_q);

G_c_R=G_I_R*camR;

G_c_T=G_I_p+G_I_R*camT;


C_p_f=G_c_R'*(G_p_f_k-G_c_T);

C_p_f=C_p_f/C_p_f(3);


Gamma=[1,0,-pts_n_temp(1);...
       0,1,-pts_n_temp(2)];


E=Gamma*C_p_f;


end










