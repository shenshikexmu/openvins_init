clear all
clc



R = [0.00181059 0.00855209 -0.99996179;
0.0414996 0.99910134 0.00861987;
0.99913688 -0.041513622 0.00145406]
q = rotMat2qRichard(R)


R_=quatern2rotMat(q);


%R=[0,1,0;1,0,0;0,0,-1]


%R=[1,0,0;0,-1,0;0,0,-1]

%R=[0,-1,0;-1,0,0;0,0,-1]

%R=[0,0,-1;...
%   0,-1,0;...
%    -1,0,0]
%
%
%%R=quatern2rotMat([-1,0,0,0])
%
%
%det(R)
%
%
%q = rotMat2qRichard_try(R);
%
%
%R_ = quatern2rotMat(q)
%
%
%
%q2=rotMat2quatern2(R);
%clc
%
%R2_=quatern2rotMat(q2)




% for i=1:1000000
%     
%     axis=randn(3,1);
%     axis=axis/norm(axis);
%     
%     q=[1;axis]*sqrt(2)*0.5;
%     
%     q=q/norm(q);
%     
%     if q(1)<0
%         q=-q;
%     end
%     
%     R=quatern2rotMat(q);
%     
%     qq=rotMat2qRichard(R);
%     
%     i
%     
%     if norm(qq-q)>1.0e-12
%         if  norm(qq-q)<(2-1.0e-12)
%             norm(qq-q)
%             
%             a=10;
%         end
%     end
%     
%   
%     
%     
%     
% end
% 
% 
% 
% 
