function imgpyr=cv_buildOpticalFlowPyramid(img,win_size,pyr_levels)

% win_size=[15,15];
% pyr_levels=5;

global matlab_or_octave

if nargin < 2
    win_size=[15,15];
end
if nargin < 3
    pyr_levels=5;
end


imgpyr{1}=img;

for i=2:pyr_levels

 %imgpyr{i,1}=cv_resize(imgpyr{i-1,1}, [0,0],0.5,0.5, 'INTER_LINEAR');
 
 [old_row,old_col,~]=size(imgpyr{i-1,1});
 
 new_row=floor(old_row*0.5);
 new_col=floor(old_col*0.5);

 if (matlab_or_octave==1)

     imgpyr{i,1}=imresize(imgpyr{i-1,1},[new_row,new_col],'bilinear');

 else
 
     imgpyr{i,1}=imresize(imgpyr{i-1,1},[new_row,new_col],"Method","linear");

 end
 

end

  
end
