function [imuPropagate]=repropagate_VINS_Mono(ImuData,Ba,Bg)
    

    
delta_p=[0;0;0];
delta_q=[1;0;0;0];    
delta_v=[0;0;0];

jacobian=eye(15);

covariance=zeros(15,15);

for i=2:size(ImuData,1)
    
    dt=ImuData(i,1)-ImuData(i-1,1);
    
    acc0=ImuData(i-1,5:7)';
    gyr0=ImuData(i-1,2:4)';
    
    acc1=ImuData(i,5:7)';
    gyr1=ImuData(i,2:4)';
   
    [delta_p,delta_q,delta_v,Ba,Bg,jacobian,covariance]=...
     midPointIntegration_VINS_Mono(dt,acc1,gyr1,acc1,gyr1,delta_p,delta_q,delta_v,Ba,Bg,jacobian,covariance);

end

sum_dt=ImuData(end,1)-ImuData(1,1);


imuPropagate{1}=delta_p;
imuPropagate{2}=delta_q;
imuPropagate{3}=delta_v;
imuPropagate{4}=Ba;
imuPropagate{5}=Bg;
imuPropagate{6}=sum_dt;
imuPropagate{7}=jacobian;
imuPropagate{8}=covariance;
    

    
end