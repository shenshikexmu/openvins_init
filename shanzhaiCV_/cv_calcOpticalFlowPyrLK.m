function [pts1, status, err] = cv_calcOpticalFlowPyrLK(imgpyr0, imgpyr1, pts0, pts1,win_size, pyr_levels, criteria, flags)

max_iters=criteria.max_iters;
eps=criteria.epsilon;

if (flags==1)
    pts1=pts1;
else
    pts1=pts0;
end

for n=pyr_levels:-1:1
    
    I0_now=imgpyr0{n};
    I1_now=imgpyr1{n};
    
    scale=2^(n-1);
    
    pts0_now=pts0/scale;
    pts1_now=pts1/scale;
    
%    figure
%    imshow([I0_now;I1_now]);
%    hold on;
%    plot(pts0_now(:,1),pts0_now(:,2),'r.');
%     
    
    [pts1] = iterOpticalflow(I0_now,I1_now,pts0_now,pts1_now,win_size,eps,max_iters);
    
    pts1=pts1*scale;
    
end

status=zeros(size(pts0,1),1);
err=zeros(size(pts0,1),1);

hsz = floor(win_size(1)/2);

[rows, cols, ~] = size(imgpyr0{1});

[x,y] = meshgrid(1:cols, 1:rows);
I0=mat2gray(imgpyr0{1});
I1=mat2gray(imgpyr1{1});

for i=1:size(pts0,1)
    
    if (pts0(i,1)<0 || pts0(i,1)>cols || pts0(i,2)<0 || pts0(i,2)>rows ||...
        pts1(i,1)<0 || pts1(i,1)>cols || pts1(i,2)<0 || pts1(i,2)>rows )
        
        err(i,1)=Inf;        
        continue;

    end

    pts0_tmp=round(pts0(i,1:2)+[0.5,0.5]);

    left = pts0_tmp(1)-hsz; right = pts0_tmp(1)+hsz;
    top = pts0_tmp(2)-hsz; bottom = pts0_tmp(2)+hsz;
    if(left<=0), left=1; end
    if(right>cols), right=cols; end
    if(top<=0), top = 1; end
    if(bottom>rows), bottom=rows; end
    win1 = I0(top:bottom,left:right);
    
    xp = x+pts1(i,1)-pts0_tmp(1);
    yp = y+pts1(i,2)-pts0_tmp(2);

    win2 = interp2(x,y,I1,xp(top:bottom,left:right),yp(top:bottom,left:right));

    it = win2-win1;
    err(i,1)=mean(mean(abs(it(~isnan(it)))));
    
    if err(i,1)<0.2
        
        status(i,1)=1;
    end

   
%    figure
%    imshow([imgpyr0{1},imgpyr1{1}]);
%    hold on;
%    plot(pts0(i,1),pts0(i,2),'r*',pts1(i,1)+size(imgpyr0{1},2),pts1(i,2),'g*');
%    hold on;
%    plot([pts0(i,1),pts1(i,1)+size(imgpyr0{1},2)],[pts0(i,2),pts1(i,2)],'b-');
%    hold off;
%    
%    
%    a=10;
    
    
end

%a=10;
                               
                               
end






function [pts1] = iterOpticalflow(I1,I2,pts0,pts1,windowSize,eps,max_iters)
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
%    change form  https://github.com/zyinshi/Optical-Flow
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if size(pts0,1)<1
    return;    
end


if size(I1,3) == 3
        I1 = rgb2gray(I1);
end
I1 = mat2gray(I1);
if size(I2,3) == 3
        I2 = rgb2gray(I2);
end
I2 = mat2gray(I2);

% h = fspecial('gaussian', [3,3], 1);
% I1 = imfilter(I1,h);
% I2 = imfilter(I2,h);

d = 1/12.*[-1,8,0,-8,1];
sz = size(I1);
Ix = conv2(I1, d,'same');
Iy = conv2(I1,d', 'same');

X2 = conv2(Ix.^2, ones(windowSize(1)),'same');
Y2 = conv2(Iy.^2, ones(windowSize(1)),'same');
XY = conv2(Ix.*Iy, ones(windowSize(1)),'same');

[x,y] = meshgrid(1:size(I1,2), 1:size(I1,1));
%% iterative

hsz = floor(windowSize(1)/2);

for i= 1:size(pts0,1)
    

    pts0_tmp=round(pts0(i,1:2)+[0.5,0.5]);

    left = pts0_tmp(1)-hsz; right = pts0_tmp(1)+hsz;
    top = pts0_tmp(2)-hsz; bottom = pts0_tmp(2)+hsz;
    if(left<=0), left=1; end
    if(right>sz(2)), right=sz(2); end
    if(top<=0), top = 1; end
    if(bottom>sz(1)), bottom=sz(1); end
    win1 = I1(top:bottom,left:right);
    ix = Ix(top:bottom,left:right);
    iy = Iy(top:bottom,left:right);
%     if i==39
%         a=10;
%     end
    A = [X2(pts0_tmp(2),pts0_tmp(1)) XY(pts0_tmp(2),pts0_tmp(1)); XY(pts0_tmp(2),pts0_tmp(1)) Y2(pts0_tmp(2),pts0_tmp(1))];
%         lambda = eig(A);
    r = rank(A);
    if(r~=2) 
        Ainv = zeros(2);
    else 
        Ainv = inv(A);
    end
    u=pts1(i,1)-pts0_tmp(1);
    v=pts1(i,2)-pts0_tmp(2);
%     mean_win1=mean(mean(win1));
    for iter = 1:max_iters
        
        xp = x+u;
        yp = y+v;

        win2 = interp2(x,y,I2,xp(top:bottom,left:right),yp(top:bottom,left:right));
%         mean_win2=mean(mean(win2));
        
%         it = win2-win1-mean_win2+mean_win1;
        it = win2-win1;
        ixt = it.*ix;
        iyt = it.*iy;
        B = -1.*[sum(ixt(:)); sum(iyt(:))];
        U = Ainv*B;
        U(isnan(U)) = 0 ;
        u = u+U(1); v = v+U(2);
        if(abs(U(1))<eps && abs(U(2))<eps) 
            break;
        end
    end
   
    pts1(i,1:2)=pts0(i,1:2)+[u,v];

%   figure;
%   
%   subplot(1,2,1)
%   imshow([I1]);
%   hold on;
%   plot(pts0(i,1),pts0(i,2),'r.');
%   title([num2str(i),'[',num2str(u),',',num2str(v),']']);
%   
%   subplot(1,2,2)
%   imshow([I2]);
%   hold on;
%   plot(pts1(i,1),pts1(i,2),'r.');
%   
%    halfwin=15;
%    for m=1:(halfwin*2-1)
%       for n=1:(halfwin*2-1)
%           
%           [x,y] = meshgrid(1:size(I1,2), 1:size(I1,1));
%            xp = x+m-halfwin;
%            yp = y+n-halfwin;
%
%            win2 = interp2(x,y,I2,xp(top:bottom,left:right),yp(top:bottom,left:right));
%       
%            it = win2-win1;
%            it_all=it.*it;
%            
%            AAA(m,n)=sum(sum(it_all));
%           
%           
%       end
%   end
%   
%   [minValue, linearIdx] = min(AAA(:));
%   [numRows, numCols] = size(AAA);
%   rowIdx = mod(linearIdx - 1, numRows) + 1;
%   colIdx = ceil(linearIdx / numRows);
%   uv=pts0(i,1:2)+[rowIdx-halfwin,colIdx-halfwin];
%   hold on;
%   plot(uv(1),uv(2),'g.');
%   title(['[',num2str(rowIdx-halfwin),',',num2str(colIdx-halfwin),']']);
%   
%
%    
%    a=10;
   
    
    
end


end     