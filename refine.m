function data2=refine(data1,status)

n=0;
data2=data1;

for i=1:size(status,1)

    if status(i,1)==1

        n=n+1;

        data2(n,:)=data1(i,:);

    end

end


data2=data2(1:n,:);


end