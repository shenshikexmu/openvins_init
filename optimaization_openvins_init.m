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


options=optimset('TolX',1e-6,'TolFun',1e-6,'Algorithm','Levenberg-Marquardt','Display','iter','MaxIter',50);

[a,resnorm]=lsqnonlin(@loss_function,a0,[],[],options,data);



for j=1:size_ids

    G_p_f_opti(:,j)=a(size_frame*15+j*3-2:size_frame*15+j*3);

end


for n=1:size_frame

   a1=a(n*15-14:n*15,1);


   q_I_k_=angleAxis2Quaternion(a1(1:3))';

   x_I_k_opti(:,n)=[q_I_k_;a1(4:15)];


end


end




function E=loss_function(a,data)

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


E=[];
precision_pic=1;


for j=1:size_ids

    G_p_f_k=a(size_frame*15+j*3-2:size_frame*15+j*3);

    for n=1:size_frame

        pts_temp=pts_opti{n}(j,1:2);

        a1=a(n*15-14:n*15,1);

        P1=a1(4:6,1);       
        if n==1
            P1=[0;0;0];
        end
        angleAxis1=a1(1:3,1);

        R1=angleAxisToRotationMatrix(angleAxis1);     %  Q1=angleAxis2Quaternion(angleAxis1)';   R1_=quatern2rotMat(Q1)
        

        Erepro=evaluate_Reprojection(P1,R1,camK,camD,camR,camT,pts_temp,G_p_f_k)*precision_pic;
    
        E=[E;Erepro];


    end


end



for i=2:size_frame
    
    %a0(i*15+1:i*15+15,1)=[angleAxistemp;Ptemp;Vtemp;Batemp;Bgtemp];
    
    a1=a(i*15-29:i*15-15,1);
    a2=a(i*15-14:i*15,1);
    
    P1=a1(4:6,1);
    if i==2
        P1=[0;0;0];
    end
    angleAxis1=a1(1:3,1);
    Q1=angleAxis2Quaternion(angleAxis1)';
    V1=a1(7:9,1);
    Ba1=a1(10:12,1);
    Bg1=a1(13:15,1);
    
    P2=a2(4:6,1);
    angleAxis2=a2(1:3,1);
    Q2=angleAxis2Quaternion(angleAxis2)';
    V2=a2(7:9,1);
    Ba2=a2(10:12,1);
    Bg2=a2(13:15,1);
    
    residuals=evaluate_VINS_Mono_gravity(P1,Q1,V1,Ba1,Bg1,P2,Q2,V2,Ba2,Bg2,imuPropagate{i},gravity); 
 
    E=[E;residuals];
end
    





end






function E=evaluate_Reprojection(G_I_p,G_I_R,camK,camD,camR,camT,pts_temp,G_p_f_k)



G_c_R=G_I_R*camR;

G_c_T=G_I_p+G_I_R*camT;


C_p_f=G_c_R'*(G_p_f_k-G_c_T);

C_p_f=C_p_f/C_p_f(3);



uv_dist = distort_cv(C_p_f(1:2), camK,camD);

E=(uv_dist-pts_temp)';


end










