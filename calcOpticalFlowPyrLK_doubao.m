function [nextPts, status, err] = calcOpticalFlowPyrLK_doubao(prevPyr, nextPyr, prevPts, winSize, maxLevel, criteria, derivLambda, flags)
    % 初始化参数
    lambda1 = 1 - derivLambda;
    lambda2 = derivLambda;
    derivKernelSize = 3;
    deriv1Scale = 0.5/4;
    deriv2Scale = 0.25/4;
    derivDepth = 'single'; % 对应 CV_32F
    halfWin = [(winSize(1)-1)*0.5 (winSize(2)-1)*0.5];

    % 断言检查
    assert(maxLevel >= 0 && winSize(1) > 2 && winSize(2) > 2);
    assert(all(size(prevPyr) == size(nextPyr)) && isequal(class(prevPyr), class(nextPyr)));

    npoints = size(prevPts, 1);
    nextPts = zeros(npoints, 2);
    status = true(npoints, 1);
    err = zeros(npoints, 1);

    if npoints == 0
        return;
    end

%     % 构建图像金字塔
%     prevPyr = cell(maxLevel+1, 1);
%     nextPyr = cell(maxLevel+1, 1);
%     prevPyr{1} = prevImg;
%     nextPyr{1} = nextImg;
%     for level = 2:maxLevel+1
%         prevPyr{level} = impyramid(prevPyr{level-1}, 'reduce');
%         nextPyr{level} = impyramid(nextPyr{level-1}, 'reduce');
%     end

    % 设置终止条件
    if ~ismember('count', criteria.type)
        criteria.maxCount = 30;
    else
        criteria.maxCount = min(max(criteria.maxCount, 0), 100);
    if ~ismember('eps', criteria.type)
        criteria.epsilon = 0.01;
    else
        criteria.epsilon = min(max(criteria.epsilon, 0), 10);
        criteria.epsilon = criteria.epsilon^2;
        
    end


    if strcmp(criteria.type , "CV_TERMCRIT_ITER")
        max_iters = criteria.max_iter;
    elseif strcmp(criteria.type , "CV_TERMCRIT_EPS")
        eps = criteria.epsilon;
    elseif strcmp(criteria.type , "CV_TERMCRIT_ITER+EPS")
        eps = criteria.epsilon;
        max_iters = criteria.max_iter;
    end








    for level = maxLevel:-1:0
        imgSize = size(prevPyr{level+1});
        % 计算导数图像缓冲区
        derivIBuf = zeros(imgSize(1)+winSize(1)*2, imgSize(2)+winSize(2)*2, 6*size(prevImg, 3), derivDepth);
        derivJBuf = zeros(imgSize(1)+winSize(1)*2, imgSize(2)+winSize(2)*2, 3*size(prevImg, 3), derivDepth);
        tempDerivBuf = zeros(imgSize, derivDepth, size(prevImg, 3));
        derivIWinBuf = zeros(winSize, derivDepth, 6*size(prevImg, 3));

        % 将当前层图像转换为指定深度并混合通道
        tempDeriv = im2single(prevPyr{level+1});
        fromTo = zeros(1, 2*size(prevImg, 3));
        for k = 1:size(prevImg, 3)
            fromTo(2*k-1) = k;
        end
        for k = 1:size(prevImg, 3)
            fromTo(2*k) = k*6;
        end
        tempDeriv = reshape(tempDeriv, [imgSize(1), imgSize(2), size(prevImg, 3)]);
        derivI = reshape(derivIBuf(winSize(1)+1:end-winSize(1), winSize(2)+1:end-winSize(2), :), [imgSize(1), imgSize(2), size(derivIBuf, 3)]);
        mixChannels(tempDeriv, derivI, fromTo);

        % 计算空间导数并合并通道
        for k = 1:size(prevImg, 3)
            tempDeriv = imfilter(prevPyr{level+1}, [1 0 -1], 'replicate', 'same', 'conv');
            fromTo(2*k) = k*6 + 1;
            mixChannels(tempDeriv, derivI, fromTo);

            tempDeriv = imfilter(prevPyr{level+1}, [1 0 -1]', 'replicate', 'same', 'conv');
            fromTo(2*k) = k*6 + 2;
            mixChannels(tempDeriv, derivI, fromTo);

            tempDeriv = imfilter(prevPyr{level+1}, [1 -2 1], 'replicate', 'same', 'conv');
            fromTo(2*k) = k*6 + 3;
            mixChannels(tempDeriv, derivI, fromTo);

            tempDeriv = imfilter(prevPyr{level+1}, [1 -1; -1 1], 'replicate', 'same', 'conv');
            fromTo(2*k) = k*6 + 4;
            mixChannels(tempDeriv, derivI, fromTo);

            tempDeriv = imfilter(prevPyr{level+1}, [1 -2 1]', 'replicate', 'same', 'conv');
            fromTo(2*k) = k*6 + 5;
            mixChannels(tempDeriv, derivI, fromTo);
        end

        tempDeriv = im2single(nextPyr{level+1});
        for k = 1:size(prevImg, 3)
            fromTo(2*k) = k*3;
        end
        mixChannels(tempDeriv, reshape(derivJBuf(winSize(1)+1:end-winSize(1), winSize(2)+1:end-winSize(2), :), [imgSize(1), imgSize(2), size(derivJBuf, 3)]), fromTo);

        for k = 1:size(prevImg, 3)
            tempDeriv = imfilter(nextPyr{level+1}, [1 0 -1], 'replicate', 'same', 'conv');
            fromTo(2*k) = k*3 + 1;
            mixChannels(tempDeriv, reshape(derivJBuf(winSize(1)+1:end-winSize(1), winSize(2)+1:end-winSize(2), :), [imgSize(1), imgSize(2), size(derivJBuf, 3)]), fromTo);

            tempDeriv = imfilter(nextPyr{level+1}, [1 0 -1]', 'replicate', 'same', 'conv');
            fromTo(2*k) = k*3 + 2;
            mixChannels(tempDeriv, reshape(derivJBuf(winSize(1)+1:end-winSize(1), winSize(2)+1:end-winSize(2), :), [imgSize(1), imgSize(2), size(derivJBuf, 3)]), fromTo);
        end

        for ptidx = 1:npoints
            prevPt = prevPts(ptidx,:) * (1/(2^level));
            nextPt = zeros(1, 2);
            if level == maxLevel
                if flags == "OPTFLOW_USE_INITIAL_FLOW"
                    nextPt = nextPts(ptidx,:) * (1/(2^level));
                else
                    nextPt = prevPt;
                end
            else
                nextPt = nextPts(ptidx,:)*2;
            end
            nextPts(ptidx,:) = nextPt;

            iprevPt = floor(prevPt) + 1; % 加一是为了和 C++代码中的索引一致
            if iprevPt(1) < -winSize(1) || iprevPt(1) >= size(derivI, 2) ||...
               iprevPt(2) < -winSize(2) || iprevPt(2) >= size(derivI, 1)
                if level == 0
                    status(ptidx) = false;
                    err(ptidx) = inf;
                end
                continue;
            end

            a = prevPt(1) - iprevPt(1);
            b = prevPt(2) - iprevPt(2);
            w00 = (1 - a)*(1 - b);
            w01 = a*(1 - b);
            w10 = (1 - a)*b;
            w11 = a*b;
            stepI = size(derivI, 3);
            stepJ = size(derivJ, 3);
            cnI = size(prevImg, 3)*6;
            cnJ = size(prevImg, 3)*3;
            A11 = 0;
            A12 = 0;
            A22 = 0;
            iA11 = 0;
            iA12 = 0;
            iA22 = 0;

            % 提取图像补丁
            for y = 1:winSize(2)
                for x = 1:winSize(1)
                    srcI = derivI(iprevPt(2)+y-1, iprevPt(1)+x-1,:);
                    dstI = derivIWinBuf(y,x,:);
                    I = srcI(1)*w00 + srcI(cnI+1)*w01 + srcI(stepI+1)*w10 + srcI(stepI+cnI+1)*w11;
                    dstI(1) = I;
                    Ix = srcI(2)*w00 + srcI(cnI+2)*w01 + srcI(stepI+2)*w10 + srcI(stepI+cnI+2)*w11;
                    Iy = srcI(3)*w00 + srcI(cnI+3)*w01 + srcI(stepI+3)*w10 + srcI(stepI+cnI+3)*w11;
                    dstI(2) = Ix; dstI(3) = Iy;
                    Ixx = srcI(4)*w00 + srcI(cnI+4)*w01 + srcI(stepI+4)*w10 + srcI(stepI+cnI+4)*w11;
                    Ixy = srcI(5)*w00 + srcI(cnI+5)*w01 + srcI(stepI+5)*w10 + srcI(stepI+cnI+5)*w11;
                    Iyy = srcI(6)*w00 + srcI(cnI+6)*w01 + srcI(stepI+6)*w10 + srcI(stepI+cnI+6)*w11;
                    dstI(4) = Ixx; dstI(5) = Ixy; dstI(6) = Iyy;

                    iA11 = iA11 + Ix*Ix;
                    iA12 = iA12 + Ix*Iy;
                    iA22 = iA22 + Iy*Iy;

                    A11 = A11 + Ixx*Ixx + Ixy*Ixy;
                    A12 = A12 + Ixy*(Ixx + Iyy);
                    A22 = A22 + Ixy*Ixy + Iyy*Iyy;
                end
            end

            A11 = lambda1*iA11 + lambda2*A11;
            A12 = lambda1*iA12 + lambda2*A12;
            A22 = lambda1*iA22 + lambda2*A22;

            D = A11*A22 - A12*A12;
            minEig = (A22 + A11 - sqrt((A11-A22)*(A11-A22) + 4*A12*A12))/(2*winSize(1)*winSize(2));
            err(ptidx) = minEig;

            if abs(D) < eps
                if level == 0
                    status(ptidx) = false;
                end
                continue;
            end

            D = 1/D;

            nextPt = nextPt - halfWin;
            prevDelta = zeros(1, 2);

            for j = 1:criteria.maxCount
                inextPt = floor(nextPt) + 1; % 加一是为了和 C++代码中的索引一致
                if inextPt(1) < -winSize(1) || inextPt(1) >= size(derivJ, 2) ||...
                   inextPt(2) < -winSize(2) || inextPt(2) >= size(derivJ, 1)
                    if level == 0
                        status(ptidx) = false;
                    end
                    break;
                end

                a = nextPt(1) - inextPt(1);
                b = nextPt(2) - inextPt(2);
                w00 = (1 - a)*(1 - b);
                w01 = a*(1 - b);
                w10 = (1 - a)*b;
                w11 = a*b;
                b1 = 0;
                b2 = 0;
                ib1 = 0;
                ib2 = 0;

                for y = 1:winSize(2)
                    for x = 1:winSize(1)
                        srcJ = derivJ(inextPt(2)+y-1, inextPt(1)+x-1,:);
                        srcI = derivIWinBuf(y,x,:);
                        It = srcJ(1)*w00 + srcJ(cnJ+1)*w01 + srcJ(stepJ+1)*w10 + srcJ(stepJ+cnJ+1)*w11 - srcI(1);
                        Ixt = srcJ(2)*w00 + srcJ(cnJ+2)*w01 + srcJ(stepJ+2)*w10 + srcJ(stepJ+cnJ+2)*w11 - srcI(2);
                        Iyt = srcJ(3)*w00 + srcJ(cnJ+3)*w01 + srcJ(stepJ+3)*w10 + srcJ(stepJ+cnJ+3)*w11 - srcI(3);
                        b1 = b1 + Ixt*srcI(4) + Iyt*srcI(5);
                        b2 = b2 + Ixt*srcI(5) + Iyt*srcI(6);
                        ib1 = ib1 + It*srcI(2);
                        ib2 = ib2 + It*srcI(3);
                    end
                end

                b1 = lambda1*ib1 + lambda2*b1;
                b2 = lambda1*ib2 + lambda2*b2;
                delta = [-((A12*b2 - A22*b1) * D), ((A12*b1 - A11*b2) * D)];
                %delta = -delta;

                nextPt = nextPt + delta;
                nextPts(ptidx,:) = nextPt + halfWin;

                if dot(delta, delta) <= criteria.epsilon
                    break;
                end

                if j > 0 && abs(delta(1)+prevDelta(1)) < 0.01 && abs(delta(2)+prevDelta(2)) < 0.01
                    nextPts(ptidx,:) = nextPts(ptidx,:) - delta*0.5;
                    break;
                end
                prevDelta = delta;
            end
        end
    end
end