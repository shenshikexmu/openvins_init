function IMU_data_fragment=get_Imu_Fragment_tc1_tc2(IMU_data,tc1,tc2)

if tc1==tc2

    IMU_data_fragment=[];

else

    for i=1:size(IMU_data,1)
    
        if IMU_data(i,1)<=tc1 && IMU_data(i+1,1)>tc1
      
           n1= i+1;
           
           a1= (IMU_data(i+1,1)-tc1)/(IMU_data(i+1,1)-IMU_data(i,1));
           
           b1= (tc1-IMU_data(i,1))/(IMU_data(i+1,1)-IMU_data(i,1));
            
        end
            
        if i>1
            
            if IMU_data(i-1,1)<tc2 && IMU_data(i,1)>=tc2
    
               n2= i-1;
    
               a2=(IMU_data(i,1)-tc2)/(IMU_data(i,1)-IMU_data(i-1,1));
    
               b2=(tc2-IMU_data(i-1,1))/(IMU_data(i,1)-IMU_data(i-1,1));
    
            end
    
        end
    
        
    end

    
    IMU_data_fragment=[IMU_data(n1-1,:)*a1+IMU_data(n1,:)*b1;...     
                       IMU_data(n1:n2,:);...
                       IMU_data(n2,:)*a2+IMU_data(n2+1,:)*b2];
    
    
                   
    % IMU_data_fragment=[IMU_data(n1,:);...     
    %                    IMU_data(n1:n2,:);...
    %                    IMU_data(n2,:)];

end

end