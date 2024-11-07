function mask=getMask(image)
  

  [rows, cols, channels] = size(image);
  
  dataType = class(image);
  
  %mask = zeros(rows, cols, 1);   
  
  mask = cast(zeros(rows,cols),dataType);


  
end
