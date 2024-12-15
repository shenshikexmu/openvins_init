function  [a,resnorm]=Optimize_my_LM2(Loss_fun,plus_fun,a0,data,TolX,TolFun,MaxIter,ConstantValue)

% author  Zhang Xin


if nargin<7
    ConstantValue=[];

end

len_uncertain=length(a0)-length(ConstantValue);



tao=1e-10;
xk=a0;

v=2;


[Ek,Jacobi]=Loss_fun(xk,data);
Jacobi=change_Jacobi_using_ConstantValue(Jacobi,ConstantValue);

g=Jacobi'*Ek;

found=logical(norm(g)<=TolFun);
mou=tao*max(diag(Jacobi'*Jacobi));

k=0;
fprintf('%12s  %12s %12s %12s \n','Iterations','Residual','Lambda','Step');
while (~found && k<MaxIter+1)
    
    %delta_x=-(Jacobi'*Jacobi+mou*eye(len_uncertain))\Jacobi'*Ek;     
    
    delta_x=-[Jacobi;sqrt(mou)*eye(len_uncertain)]\[Ek;zeros(len_uncertain,1)];  

    %     M=-(Jacobi'*Jacobi+Lambda*0.001*sqrt(diag(diag(Jacobi'*Jacobi)))*eye(len_uncertain));
%     P=Jacobi'*Ek;
%     delta_x=schur_complement(M,P,data);

    
    if (norm(delta_x)<=TolX*(norm(xk)+TolX))
        found=true;

    else
        xk_new=xk_plus_delta_x(plus_fun,xk,delta_x,ConstantValue,data);%xk_new=xk+delta_x';
%         Ek=Loss_fun(xk,data);
%         Ek_new=Loss_fun(xk_new,data);

        [Ek_new,Jacobi_new]=Loss_fun(xk_new,data);
        Jacobi_new=change_Jacobi_using_ConstantValue(Jacobi_new,ConstantValue);

        L0=delta_x'*Jacobi'*Ek;
        L_delta=delta_x'*Jacobi'*Jacobi*delta_x;
        rho=(Ek'*Ek-Ek_new'*Ek_new)/(-L0-L_delta);
        
        
        if rho>0

            fprintf('%7d  %18f %12f %15.8f \n',k,Ek'*Ek, mou, norm(delta_x));
            
            %fprintf('Iterations: %d, Residual: %d, Step: %d \n',k,round(Ek'*Ek,2),norm(delta_x));
            k=k+1;
            
            xk=xk_new;
            %Jacobi=Get_Jacobi(Loss_fun,xk,data);
            %Ek=Loss_fun(xk,data);
            Jacobi=Jacobi_new;
            Ek=Ek_new;
            
            g=Jacobi'*Ek;
            found=(norm(g)<=TolFun);
            mou=mou*max([1/3,1-(2*rho-1)^3]);
            v=2;
        else
            mou=mou*v;
            v=2*v;
           
        end
  
    end

end


% xk=xk_plus_delta_x(plus_fun,xk,delta_x,ConstantValue,data); %xk=xk+delta_x';
% Ek=Loss_fun(xk,data);
% fprintf('%7d  %18f %12f %15.8f \n',k, Ek'*Ek, mou, norm(delta_x));



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


function Jacobi=Get_Jacobi(Loss_fun,xk,data)

scale=1e-4;

%Ek=Loss_fun(xk,data);

for i=1:length(xk)
    x_temp1=xk;
    x_temp2=xk;
    if abs(x_temp1(i))>scale
  
        delta=x_temp1(i)*scale;
    else
        delta=scale;
    end
    x_temp1(i)=x_temp1(i)+delta;
    x_temp2(i)=x_temp2(i)-delta;

    E_temp1=Loss_fun(x_temp1,data);

    E_temp2=Loss_fun(x_temp2,data);

    Jacobi(:,i)=(E_temp1-E_temp2)/delta/2;

end

end