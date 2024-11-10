function [nextPts, status, err] = calcOpticalFlowPyrLK(prevImg, nextImg, prevPts, nextPts, winSize, maxLevel, criteria, derivLambda, flags)

    derivLambda = min(max(derivLambda, 0), 1);
    lambda1 = 1 - derivLambda;
    lambda2 = derivLambda;
    derivKernelSize = 3;
    deriv1Scale = 0.5 / 4.0;
    deriv2Scale = 0.25 / 4.0;
    derivDepth = 'single';  % In MATLAB, we use 'single' for float
    halfWin = [(winSize(1)-1)*0.5, (winSize(2)-1)*0.5];
    
    % Ensure proper sizes
    assert(maxLevel >= 0 && winSize(1) > 2 && winSize(2) > 2);
    assert(all(size(prevImg) == size(nextImg)));
    
    % Initialize points
    npoints = size(prevPts, 1);
    nextPts = zeros(npoints, 2);
    status = true(npoints, 1);
    err = zeros(npoints, 1);

    if npoints == 0
        return;
    end
    
    % Build image pyramids
    prevPyr = cell(maxLevel + 1, 1);
    nextPyr = cell(maxLevel + 1, 1);
    prevPyr{1} = prevImg;
    nextPyr{1} = nextImg;
    
    for i = 2:maxLevel + 1
        prevPyr{i} = imresize(prevPyr{i-1}, 0.5);
        nextPyr{i} = imresize(nextPyr{i-1}, 0.5);
    end
    
    % Initialize derivatives (MATLAB does not have the same type system as OpenCV)
    for level = maxLevel:-1:0
        imgSize = size(prevPyr{level+1});
        [Ix, Iy] = imgradientxy(prevPyr{level+1}, 'sobel');
        
        for ptidx = 1:npoints
            prevPt = prevPts(ptidx, :) / (2^level);
            
            if level == maxLevel
                if bitand(flags, 1)  % OPTFLOW_USE_INITIAL_FLOW
                    nextPt = nextPts(ptidx, :) / (2^level);
                else
                    nextPt = prevPt;
                end
            else
                nextPt = nextPts(ptidx, :) * 2;
            end
            
            % Check boundaries
            prevPt = prevPt - halfWin;
            iprevPt = floor(prevPt);
            
            if iprevPt(1) < 1 || iprevPt(1) >= size(Ix, 2) || iprevPt(2) < 1 || iprevPt(2) >= size(Iy, 1)
                if level == 0
                    status(ptidx) = false;
                    err(ptidx) = inf;
                end
                continue;
            end
            
            % Iterate over levels and apply the Lucas-Kanade iterations
            for j = 1:criteria.maxCount
                nextPt = nextPt - halfWin;
                inextPt = floor(nextPt);
                
                if inextPt(1) < 1 || inextPt(1) >= size(Ix, 2) || inextPt(2) < 1 || inextPt(2) >= size(Iy, 1)
                    if level == 0
                        status(ptidx) = false;
                    end
                    break;
                end
                
                % Compute the update for the flow vector based on image gradients
                delta = ...  % compute flow update here based on A matrix and gradients

                % Compute flow update delta based on A and b
                detA = A11 * A22 - A12^2;  % Determinant of A
                if abs(detA) < eps
                    delta = [0; 0];  % If determinant is too small, skip update
                else
                    invA = [A22, -A12; -A12, A11] / detA;  % Inverse of A
                    delta = -invA * [b1; b2];  % Compute delta (flow update)
                end
                
                nextPt = nextPt + delta;
                nextPts(ptidx, :) = nextPt + halfWin;
                
                if sum(delta.^2) <= criteria.epsilon^2
                    break;
                end
            end
        end
    end
end
