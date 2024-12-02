function  [a,resnorm]=Optimize_my_LM(Loss_fun,a0,data,TolX,TolFun,MaxIter,ConstantValue)

% author  Zhang Xin

if nargin<7
    ConstantValue=[];

end

len_uncertain=length(a0)-length(ConstantValue);


Lambda=1e-2;
xk=a0;




Jacobi=Get_Jacobi(Loss_fun,xk,data,ConstantValue);
Ek=Loss_fun(xk,data);

g=Jacobi'*Ek;

found=logical(norm(g)<=TolFun);


k=0;
fprintf('%12s  %12s %12s %12s \n','Iterations','Residual','Lambda','Step');
while (~found && k<MaxIter+1)
    
    delta_x=-(Jacobi'*Jacobi+Lambda*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain))\Jacobi'*Ek;    
    
    %delta_x=-[Jacobi;(Lambda)^2*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain)]\[Ek;zeros(len_uncertain,1)];      

    
    if (norm(delta_x)<=TolX*(norm(xk)+TolX))
        found=true;

    else
        xk_new=xk_plus_delta_x(xk,delta_x,ConstantValue);%xk_new=xk+delta_x';
        Ek=Loss_fun(xk,data);
        Ek_new=Loss_fun(xk_new,data);
        L0=delta_x'*Jacobi'*Ek;
        L_delta=delta_x'*Jacobi'*Jacobi*delta_x;
        rho=(Ek'*Ek-Ek_new'*Ek_new)/(-L0-L_delta);
        
        
        if rho>0
            
            fprintf('%7d  %18d %12f %15.8f \n',k, Ek'*Ek, Lambda, norm(delta_x));
            k=k+1;
            
            found=(norm(Ek'*Ek-Ek_new'*Ek_new)<=TolFun);  
            xk=xk_new;
            Jacobi=Get_Jacobi(Loss_fun,xk,data,ConstantValue);
            Ek=Loss_fun(xk,data);
       
            Lambda=Lambda/10;
        
        else
            Lambda=Lambda*10;
            
        end
  
    end

end


xk=xk_plus_delta_x(xk,delta_x,ConstantValue); %xk=xk+delta_x';
Ek=Loss_fun(xk,data);
fprintf('%7d  %18d %12f %15.8f \n',k, Ek'*Ek, Lambda, norm(delta_x));


a=xk;

resnorm=Ek'*Ek;


end



function xk=xk_plus_delta_x(xk,delta_x,ConstantValue)


n_ConstantValue=0;

for i=1:size(xk)


    if i_in_ConstantValue(i,ConstantValue)

        n_ConstantValue=n_ConstantValue+1;

    else
        xk(i)=xk(i)+delta_x(i-n_ConstantValue);
    end

end



end


function Jacobi=Get_Jacobi(Loss_fun,xk,data,ConstantValue)

scale=1e-4;

Ek=Loss_fun(xk,data);


n_ConstantValue=0;

for i=1:length(xk)
    
    if i_in_ConstantValue(i,ConstantValue)
        n_ConstantValue=n_ConstantValue+1;
        continue;
    end

    x_temp1=xk;
   % x_temp2=xk;
    if abs(x_temp1(i))>scale
  
        delta=x_temp1(i)*scale;
    else
        delta=scale;
    end
    x_temp1(i)=x_temp1(i)+delta;
  %  x_temp2(i)=x_temp2(i)-delta;

    E_temp1=Loss_fun(x_temp1,data);

  %  E_temp2=Loss_fun(x_temp2,data);

    Jacobi(:,i-n_ConstantValue)=(E_temp1-Ek)/delta;

end

end





function bool=i_in_ConstantValue(i,ConstantValue)

bool=0;

if isempty(ConstantValue)

    return ;


else

    for j=1:length(ConstantValue)

        if i==ConstantValue(j)

            bool=1;
            break;
        end

    end

end

end