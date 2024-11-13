function img_out=cv_equalizeHist(img_in)


global matlab_or_octave


if (matlab_or_octave==1)      % matlab

    img_out=histeq(img_in);

else                          % octave
  
    img_in=histeq(img_in);
    img_out=uint8(img_in*255);

end
    

end
