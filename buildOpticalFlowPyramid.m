function imgpyr=buildOpticalFlowPyramid(img,win_size,pyr_levels)

  % win_size=[15,15];
  % pyr_levels=5;
  
  if nargin < 2
      win_size=[15,15];
  end
  if nargin < 3
      pyr_levels=5;
  end


  imgpyr{1}=img;

  scaleFactor=0.5;

  for i=2:pyr_levels
    
     imgpyr{i}=cv_resize(imgpyr{i-1}, [0,0],0.5,0.5, 'INTER_LINEAR');
    
  end

  
end
