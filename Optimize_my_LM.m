function  [a,resnorm]=Optimize_my_LM(Loss_fun,plus_fun,a0,data,TolX,TolFun,MaxIter,ConstantValue)

% author  Zhang Xin

if nargin<7
    ConstantValue=[];

end

len_uncertain=length(a0)-length(ConstantValue);




Lambda=1e-2;
xk=a0;


[Ek,Jacobi]=Loss_fun(xk,data);
Jacobi=change_Jacobi_using_ConstantValue(Jacobi,ConstantValue);

%Jacobibi=Get_Jacobi(Loss_fun,plus_fun,xk,data,ConstantValue);



g=Jacobi'*Ek;

found=logical(norm(g)<=TolFun);


k=0;
fprintf('%12s  %12s %12s %12s \n','Iterations','Residual','Lambda','Step');
while (~found && k<MaxIter+1)
    
    delta_x=-(Jacobi'*Jacobi+Lambda*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain))\Jacobi'*Ek;    

%     M=-(Jacobi'*Jacobi+Lambda*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain));
%     P=Jacobi'*Ek;
%     delta_x=schur_complement(M,P,data);


    
    %delta_x=-[Jacobi;(Lambda)*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain)]\[Ek;zeros(len_uncertain,1)]; 

    %delta_x=-inv(Jacobi'*Jacobi+Lambda*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain))*Jacobi'*Ek;

    
    if (norm(delta_x)<=TolX*(norm(xk)+TolX))
        found=true;

    else
        xk_new=xk_plus_delta_x(plus_fun,xk,delta_x,ConstantValue,data);%xk_new=xk+delta_x';
        
        %Ek_new=Loss_fun(xk_new,data);

        [Ek_new,Jacobi_new]=Loss_fun(xk_new,data);
        Jacobi_new=change_Jacobi_using_ConstantValue(Jacobi_new,ConstantValue);


        L0=delta_x'*Jacobi'*Ek;
        L_delta=delta_x'*Jacobi'*Jacobi*delta_x;
        rho=(Ek'*Ek-Ek_new'*Ek_new)/(-L0-L_delta);

        rho=(Ek'*Ek-Ek_new'*Ek_new);
        
        
        if rho>0
            
            fprintf('%7d  %18d %12f %15.8f \n',k, Ek'*Ek, Lambda, norm(delta_x));
            k=k+1;
            
            found=(norm(Ek'*Ek-Ek_new'*Ek_new)<=TolFun);  
            xk=xk_new;
            Jacobi=Jacobi_new;
            Ek=Ek_new;
       
            Lambda=Lambda/10;
        
        else
            Lambda=Lambda*10;
            
        end
  
    end

end


% xk=xk_plus_delta_x(plus_fun,xk,delta_x,ConstantValue,data); %xk=xk+delta_x';
% Ek=Loss_fun(xk,data);
% fprintf('%7d  %18d %12f %15.8f \n',k, Ek'*Ek, Lambda, norm(delta_x));


a=xk;

resnorm=Ek'*Ek;


end




function  X=schur_complement(M,P,data)


% M*X=P
%size_frame=data{8};
size_ids=data{9};


M_n=size(M,1);

A_n=M_n-size_ids*3;

%D_n=size_ids*3;


A=M(1:A_n,1:A_n);

a=P(1:A_n);

b=P(A_n+1:M_n);

B=M(1:A_n,A_n+1:M_n);

C=M(A_n+1:M_n,1:A_n);

D=M(A_n+1:M_n,A_n+1:M_n);

invD=inv(D);


x=(A-B*invD*C)\(a-B*invD*b);

y=B\(B*invD*b-B*invD*C*x);


X=[x;y];




end




function Jacobi_new=change_Jacobi_using_ConstantValue(Jacobi,ConstantValue)

if isempty(ConstantValue)

    Jacobi_new=Jacobi;

else

    Jacobi_new=zeros(size(Jacobi,1),size(Jacobi,2)-length(ConstantValue));

    n_ConstantValue=0;

    for i=1:size(Jacobi,2)

        if i_in_ConstantValue(i,ConstantValue)

            n_ConstantValue=n_ConstantValue+1;

        else

            Jacobi_new(:,i-n_ConstantValue)= Jacobi(:,i);

        end

    end

end


end






function xk=xk_plus_delta_x(plus_fun,xk,delta_x,ConstantValue,data)


n_ConstantValue=0;

delta_x2=zeros(size(xk));

for i=1:size(xk,1)


    if i_in_ConstantValue(i,ConstantValue)

        n_ConstantValue=n_ConstantValue+1;

    else
        delta_x2(i,1)=delta_x(i-n_ConstantValue);
    end

end


xk=plus_fun(xk,delta_x2,data);


end




function Jacobi=Get_Jacobi(Loss_fun,plus_fun,xk,data,ConstantValue)

scale=1e-3;

Ek=Loss_fun(xk,data);


n_ConstantValue=0;

for i=1:length(xk)
    
    if i_in_ConstantValue(i,ConstantValue)
        n_ConstantValue=n_ConstantValue+1;
        continue;
    end

    x_delta=xk*0;
   % x_temp2=xk;
    if abs(xk(i))>scale
  
        delta=xk(i)*scale;
    else
        delta=scale;
    end
    x_delta(i)=x_delta(i)+delta;
  %  x_temp2(i)=x_temp2(i)-delta;

    x_temp1=plus_fun(xk,x_delta,data);

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