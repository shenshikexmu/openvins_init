clear all
clc
addpath('shanzhaiCV');
addpath('yamlMatlab');
addpath('repropagate');
addpath('quaternion');

global matlab_or_octave
matlab_or_octave=1; 

    
file_cam0='../bag/V1_02_medium/mav0/cam0/';
file_imu0='../bag/V1_02_medium/mav0/imu0/';

if  (matlab_or_octave ==1)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  matlab   
    datacsv_cam0=readcell([file_cam0,'data.csv']);
    datacsv_imu0=readcell([file_imu0,'data.csv']);
    imuData=cell2mat(datacsv_imu0(2:end,:));
    imuData(:,1)=imuData(:,1)*10e-10;
    cam0Para = ReadYaml_matlab([file_cam0,'sensor.yaml']);            % cell2mat
    cam0Para.intrinsics=cell2mat(cam0Para.intrinsics);
    cam0Para.distortion_coefficients=cell2mat(cam0Para.distortion_coefficients);
    cam0Para.T_BS.data=reshape(cell2mat(cam0Para.T_BS.data),4,4)';
else                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  octave
    pkg load io
    pkg load image
    datacsv_cam0 = csv2cell([file_cam0,'data.csv']);
    cam0Para = readYaml_octave([file_cam0,'sensor.yaml']);
end


global min_px_dist  grid_x grid_y num_features threshold currid

gravity_mag= 9.81;

min_px_dist=10;
grid_x=5; 
grid_y=5;
currid=0;


pyr_levels=5;

num_pts=200;  % 定位需要提取的特征点数
num_cameras=1;

num_features=round(num_pts/num_cameras);

fast_threshold=20;

threshold=fast_threshold;

win_size = [15, 15];

camK=[cam0Para.intrinsics(1),0,cam0Para.intrinsics(3);...
      0,cam0Para.intrinsics(2),cam0Para.intrinsics(4);...
      0,0,1];

camD=[cam0Para.distortion_coefficients(1),cam0Para.distortion_coefficients(2),cam0Para.distortion_coefficients(3),cam0Para.distortion_coefficients(4)];

camR=cam0Para.T_BS.data(1:3,1:3);

camT=cam0Para.T_BS.data(1:3,4);


max_focallength=max(camK(1,1),camK(2,2));


global ACC_N GYR_N ACC_W GYR_W 


ACC_N=0.5;                     %bao
GYR_N=0.1;
ACC_W=0.001;
GYR_W=0.00001;




i=989;



m=10;
for n=1:m
    
    img{n}=imread([file_cam0,'data/',datacsv_cam0{i+n-1,2}]);
    img{n}=cv_equalizeHist(img{n});
    imgpyr{n}=cv_buildOpticalFlowPyramid(img{n},win_size,pyr_levels);
    
    if (exist(['pts_ids_',num2str(i),'_',num2str(m),'.mat'])==2)   
        if n==m           
            load(['pts_ids_',num2str(i),'_',num2str(m),'.mat']);
        end 
        continue;   
    end
      
    if n==1
        mask=getMask(img{n});
        pts{n}=[];
        ids{n}=[];
        [pts{n},ids{n}]=perform_detection_monocular(imgpyr{n},mask,pts{n},ids{n});
        continue;   
    end
    
    criteria.max_iters=30;
    criteria.epsilon=0.01;
    cv_OPTFLOW_USE_INITIAL_FLOW=1;
    
    pts{n}=pts{n-1};
    [pts{n}, status, err] = cv_calcOpticalFlowPyrLK(imgpyr{n-1}, imgpyr{n}, pts{n-1}, pts{n}, win_size, pyr_levels, criteria, cv_OPTFLOW_USE_INITIAL_FLOW);

    pts{n-1}=refine(pts{n-1},status);
    ids{n-1}=refine(ids{n-1},status);
    pts{n}=refine(pts{n},status);

    pts_n{n-1}=zeros(size(pts{n-1},1),2);
    pts_n{n}=zeros(size(pts{n},1),2);

    for p=1:size(pts{n-1},1)
        pts_n{n-1}(p,:)=undistort_cv(pts{n-1}(p,1:2)-[1,1], camK,camD);
        pts_n{n}(p,:)=undistort_cv(pts{n}(p,1:2)-[1,1], camK,camD);
    end

    [mask]=cv_findFundamentalMat(pts_n{n-1}, pts_n{n}, 'cv_FM_RANSAC', 1/max_focallength ,0.999);

    pts{n-1}=refine(pts{n-1},mask);
    pts_n{n-1}=refine(pts_n{n-1},mask);
    ids{n-1}=refine(ids{n-1},mask);
    pts{n}=refine(pts{n},mask);
    pts_n{n}=refine(pts_n{n},mask);
    ids{n}=ids{n-1};
    
    if n==m
        save(['pts_ids_',num2str(i),'_',num2str(n),'.mat'],'pts','ids','pts_n');
    end
    
end

%%

ba=[0;0;0];
bg=[0;0;0];

for n=1:m

    camTimeStamp(n,1)=datacsv_cam0{i+n-1,1}*10e-10;

    if n>1

        imuData_fragment{n}=get_Imu_Fragment_tc1_tc2(imuData,camTimeStamp(n-1),camTimeStamp(n));

    else

        imuData_fragment{n}=[];

    end

    imuPropagate{n}=repropagate_VINS_Mono(imuData_fragment{n},ba,bg);

    imuData_fragment_joint{n}=get_Imu_Fragment_tc1_tc2(imuData,camTimeStamp(1),camTimeStamp(n));

    imuPropagate_joint{n}=repropagate_VINS_Mono(imuData_fragment_joint{n},ba,bg);

    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  paper <OpenVINS State Initiaolization: Details and Derivations> 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


size_ids=size(ids{10},1);
size_frame=10;
M0=zeros(size_ids*size_frame*2,size_ids*3+6);
b0=zeros(size_ids*size_frame*2,1);

for i=1:size_ids%size(ids{10},1)

    for n=1:size_frame

        for s=1:size(ids{n},1)

             if ids{n}(s)==ids{10}(i)

                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 pts_opti{n}(i,:)=pts{n}(s,:);
                 ids_opti{n}(i,:)=ids{n}(s,:);
                 pts_n_opti{n}(i,:)=pts_n{n}(s,:);


                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 Gamma=[1,0,-pts_n{n}(s,1);...
                        0,1,-pts_n{n}(s,2)];

                 C_I_R=camR';

                 I0_alpha_k=imuPropagate_joint{n}{1};

                 I0_q_k=imuPropagate_joint{n}{2};

                 Delta_T=imuPropagate_joint{n}{6};

                 Ik_I0_R=quatern2rotMat(I0_q_k)';

                 % 论文中公式（28）（29）（30）存在错误，这里计算的时候做了修改
                 M((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,(i-1)*3+1:(i-1)*3+3)=Gamma*C_I_R*Ik_I0_R;

                 M((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,size_ids*3+1:size_ids*3+3)=-Gamma*C_I_R*Ik_I0_R*Delta_T;

                 M((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,size_ids*3+4:size_ids*3+6)=0.5*Gamma*C_I_R*Ik_I0_R*Delta_T*Delta_T;

                 b((i-1)*size_frame*2+n*2-1:(i-1)*size_frame*2+n*2,1)=Gamma*camR'*camT+ Gamma*C_I_R*Ik_I0_R*I0_alpha_k;

             end

        end

    end

end

%%
A1=M(:,1:size_ids*3+3);
A2=M(:,size_ids*3+4:size_ids*3+6);

inv_A1A1=inv(A1'*A1);

tempMatrix=eye(size_ids*size_frame*2)-A1*inv_A1A1*A1';
D=A2'*tempMatrix*A2;

d=A2'*tempMatrix*b;

lambda=getLambda(D,d,gravity_mag);

I0_g=inv(D-lambda*eye(3))*d;

x1=-inv_A1A1*A1'*A2*I0_g+inv_A1A1*A1'*b;


x=M\b;

%%


I0_v_I0=x1(size_ids*3+1:size_ids*3+3,1);

G_I0_q=getInitQuaternionfromAcc_yfront(I0_g);
G_I0_R=quatern2rotMat(G_I0_q);

for i=1:size_ids

    I0_p_f(:,i)=x1(i*3-2:i*3);

    G_p_f(:,i)=G_I0_R*I0_p_f(:,i);

end




for n=1:size_frame

    I0_Ik_q=imuPropagate_joint{n}{2};
    I0_Ik_R=quatern2rotMat(I0_Ik_q);
    delta_T=imuPropagate_joint{n}{6};
    I0_alpha_k=imuPropagate_joint{n}{1};
    I0_bata_k=imuPropagate_joint{n}{3};

    I0_p_Ik=I0_v_I0*delta_T-0.5*I0_g*delta_T^2+I0_alpha_k;
    I0_v_Ik=I0_v_I0-I0_g*delta_T+I0_bata_k;

    G_Ik_q=quaternProd(G_I0_q,I0_Ik_q);  %   G_Ik_R=quatern2rotMat(G_Ik_q)    G_Ik_R_=G_I0_R*I0_Ik_R
    G_p_Ik=G_I0_R*I0_p_Ik;
    G_v_Ik=G_I0_R*I0_v_Ik;


     x_I_k(:,n)=[G_Ik_q;G_p_Ik;G_v_Ik;bg;ba];


end


%%




[x_I_k_opti,G_p_f_opti]=optimaization_openvins_init(x_I_k,G_p_f,imuPropagate,pts_opti,pts_n_opti,camK,camD,camR,camT,gravity_mag);




%%


figure

for i=1:size(G_p_f_opti,2)

   G_p_f_k=G_p_f_opti(:,i);
   plot3(G_p_f_k(1),G_p_f_k(2),G_p_f_k(3),'r*');
   hold on;


end

img_x_len=0.2;
img_y_len=0.127;

for i=1:size(x_I_k_opti,2)


    x_I_k_=x_I_k_opti(:,i);

    G_p_k=x_I_k_(5:7);

    G_q_k=x_I_k_(1:4);

    G_R_k=quatern2rotMat(G_q_k);

    G_camR_k=G_R_k*camR;

    G_camT_k=G_p_k+G_R_k*camT;

    p1=G_camT_k-G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;
    p2=G_camT_k-G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
    p3=G_camT_k+G_camR_k(:,1)*img_x_len+G_camR_k(:,2)*img_y_len;
    p4=G_camT_k+G_camR_k(:,1)*img_x_len-G_camR_k(:,2)*img_y_len;

    lines=[p1,p2,p3,p4,p1];

    plot3(lines(1,:),lines(2,:),lines(3,:),'b-');
    hold on;

end


axis equal;

hold off;










%%

if 1

    for n=size_frame
          
        colors = jet(6);
        
        figure;
        imshow([imgpyr{1}{1},imgpyr{n}{1};imgpyr{n}{1},imgpyr{n}{1}]);
        hold on;
        
        s=0;
        for k=1:size(ids_opti{1},1)   
           
                    
            s=s+1;
            hold on
            plot(pts_opti{1}(k,1),pts_opti{1}(k,2),'r+');
            
            if mod(s,2)==0
                hold on;
                plot(pts_opti{n}(k,1),pts_opti{n}(k,2)+size(imgpyr{1}{1},1),'g+');
                hold on;
                plot([pts_opti{1}(k,1),pts_opti{n}(k,1)],[pts_opti{1}(k,2),pts_opti{n}(k,2)+size(imgpyr{1}{1},1)],'-', 'Color', colors(mod(s,6)+1,:));
            else
                hold on;
                plot(pts_opti{n}(k,1)+size(imgpyr{1}{1},2),pts_opti{n}(k,2),'g+');
                hold on;
                plot([pts_opti{1}(k,1),pts_opti{n}(k,1)+size(imgpyr{1}{1},2)],[pts_opti{1}(k,2),pts_opti{n}(k,2)],'-', 'Color', colors(mod(s,6)+1,:));
            end

            hold on;
            %quiver(pts{1}(j,1)+size(imgpyr{1}{1},2), pts{1}(j,2)+size(imgpyr{1}{1},1), pts{n}(k,1)-pts{1}(j,1), pts{n}(k,2)-pts{1}(j,2), 'r', 'LineWidth', 2);

            plot([pts_opti{1}(k,1)+size(imgpyr{1}{1},2),pts_opti{n}(k,1)+size(imgpyr{1}{1},2)],[pts_opti{1}(k,2)+size(imgpyr{1}{1},1),pts_opti{n}(k,2)+size(imgpyr{1}{1},1)], 'r', 'LineWidth', 2);
    
               
            
        end
        
        title(['Optical Flow (Lucas-Kanade) image',num2str(1),' to image',num2str(n)]);
        
        hold off;
        
        a=10;
    
        
    end

end










