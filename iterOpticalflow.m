function [pts1] = iterOpticalflow(I1,I2,pts0,pts1,windowSize,eps,max_iters)
    
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


%% iterative

hsz = floor(windowSize(1)/2);

for i= 1:size(pts0,1)
    
    pts0_tmp=round(pts0(i,:));

    left = pts0_tmp(1)-hsz; right = pts0_tmp(1)+hsz;
    top = pts0_tmp(2)-hsz; bottom = pts0_tmp(2)+hsz;
    if(left<=0), left=1; end
    if(right>sz(2)), right=sz(2); end
    if(top<=0), top = 1; end
    if(bottom>sz(1)), bottom=sz(1); end
    win1 = I1(top:bottom,left:right);
    ix = Ix(top:bottom,left:right);
    iy = Iy(top:bottom,left:right);
    if i==84
        a=10;
    end
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
    for iter = 1:max_iters
        [x,y] = meshgrid(1:size(I1,2), 1:size(I1,1));
        xp = x+u;
        yp = y+v;

        win2 = interp2(x,y,I2,xp(top:bottom,left:right),yp(top:bottom,left:right));
        
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
    i
    pts1(i,1:2)=pts0(i,1:2)+[u,v];
   [u,v]
   
   figure;
   subplot(1,2,1)
   imshow([I1]);
   hold on;
   plot(pts0(i,1),pts0(i,2),'r.');
   subplot(1,2,2)
   imshow([I2]);
   hold on;
   plot(pts1(i,1),pts1(i,2),'r.');
   
   
    for m=1:31
       for n=1:31
           
           [x,y] = meshgrid(1:size(I1,2), 1:size(I1,1));
            xp = x+m-16;
            yp = y+n-16;

            win2 = interp2(x,y,I2,xp(top:bottom,left:right),yp(top:bottom,left:right));
       
            it = win2-win1;
            it_all=it.*it;
            
            AAA(m,n)=sum(sum(it_all));
           
           
       end
   end
   
   [minValue, linearIdx] = min(AAA(:));
   [numRows, numCols] = size(AAA);
   rowIdx = mod(linearIdx - 1, numRows) + 1;
   colIdx = ceil(linearIdx / numRows);
   [rowIdx-16,colIdx-16]
   uv=pts0(i,1:2)+[rowIdx-16,colIdx-16];
   hold on;
   plot(uv(1),uv(2),'g.');
   

    
    a=10;
   
    
    
end


end     