function  [result_delta_p,result_delta_q,result_delta_v,result_linearized_ba,result_linearized_bg,result_jacobian,result_covariance]=...
          midPointIntegration_VINS_Mono(dt,acc0,gyr0,acc1,gyr1,delta_p,delta_q,delta_v,linearized_ba,linearized_bg,jacobian,covariance)
    
global ACC_N GYR_N ACC_W GYR_W ;
    
unacc0=quatern2rotMat(delta_q)*(acc0-linearized_ba);
ungyr=0.5*(gyr0+gyr1)-linearized_bg;
    
result_delta_q=[ 1           , -ungyr(1)*dt/2 , -ungyr(2)*dt/2 , -ungyr(3)*dt/2  ;...
                ungyr(1)*dt/2 ,   1           ,  ungyr(3)*dt/2 , -ungyr(2)*dt/2  ;...
                ungyr(2)*dt/2 , -ungyr(3)*dt/2 ,   1           ,  ungyr(1)*dt/2  ;...
                ungyr(3)*dt/2 ,  ungyr(2)*dt/2 , -ungyr(1)*dt/2 ,           1   ]*delta_q;   
                        
result_delta_q=result_delta_q/norm(result_delta_q);
    
unacc1=quatern2rotMat(result_delta_q)*(acc1-linearized_ba);
    
unacc=0.5*(unacc0+unacc1);
    
result_delta_p=delta_p + delta_v *dt + 0.5*unacc * dt *dt;

result_delta_v=delta_v+ unacc*dt;

result_linearized_ba = linearized_ba;
result_linearized_bg = linearized_bg;
    
wx=0.5*(gyr0+gyr1)-linearized_bg;

a0x=acc0-linearized_ba;
a1x=acc1-linearized_ba;
    
Rwx=Skew_symmetric(wx);
Ra0x=Skew_symmetric(a0x);
Ra1x=Skew_symmetric(a1x);  
    
F=zeros(15,15);
F(1:3,1:3)=eye(3);
F(1:3,4:6)=-0.25*quatern2rotMat(delta_q)*Ra0x*dt*dt-0.25*quatern2rotMat(result_delta_q)*Ra1x*(eye(3)-Rwx*dt)*dt*dt;
F(1:3,7:9)=eye(3)*dt;
F(1:3,10:12)=-0.25*(quatern2rotMat(delta_q)+quatern2rotMat(result_delta_q))*dt*dt;
F(1:3,13:15)=0.25*quatern2rotMat(result_delta_q)*Ra1x*dt*dt*dt;
F(4:6,4:6)=eye(3)-Rwx*dt;
F(4:6,13:15)=-1.0*eye(3)*dt;
F(7:9,4:6)=-0.5*quatern2rotMat(delta_q)*Ra0x*dt-0.5*quatern2rotMat(result_delta_q)*Ra1x*(eye(3)-Rwx*dt)*dt;
F(7:9,7:9)=eye(3);
F(7:9,10:12)=-0.5*(quatern2rotMat(delta_q)+quatern2rotMat(result_delta_q))*dt;
F(7:9,13:15)=0.5*quatern2rotMat(result_delta_q)*Ra1x*dt*dt;
F(10:12,10:12)=eye(3);
F(13:15,13:15)=eye(3);

V=zeros(15,18);                                              % Need to confirm, the code is not the same as the article
V(1:3,1:3)=-0.25*quatern2rotMat(delta_q)*dt*dt;
V(1:3,4:6)=0.25*quatern2rotMat(result_delta_q)*Ra1x*dt*dt*0.5*dt;
V(1:3,7:9)=-0.25*quatern2rotMat(result_delta_q)*dt*dt;
V(1:3,10:12)=0.25*quatern2rotMat(result_delta_q)*Ra1x*dt*dt*0.5*dt;
V(4:6,4:6)=0.5*eye(3)*dt;
V(4:6,10:12)=0.5*eye(3)*dt;
V(7:9,1:3)=0.5*quatern2rotMat(delta_q)*dt;
V(7:9,4:6)=-0.5*quatern2rotMat(result_delta_q)*Ra1x*dt*0.5*dt;
V(7:9,7:9)=0.5*quatern2rotMat(result_delta_q)*dt;
V(7:9,10:12)=-0.5*quatern2rotMat(result_delta_q)*Ra1x*dt*0.5*dt;
V(10:12,13:15)=eye(3)*dt;
V(13:15,16:18)=eye(3)*dt;


result_jacobian=F*jacobian;

noise=zeros(18,18);

noise(1:3,1:3)=ACC_N*ACC_N*eye(3);
noise(4:6,4:6)=GYR_N*GYR_N*eye(3);
noise(7:9,7:9)=ACC_N*ACC_N*eye(3);
noise(10:12,10:12)=GYR_N*GYR_N*eye(3);
noise(13:15,13:15)=ACC_W*ACC_W*eye(3);
noise(16:18,16:18)=GYR_W*GYR_W*eye(3);



result_covariance=F*covariance*F'+V*noise*V';
    

    
end


