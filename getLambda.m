function result=getLambda(D,d,g)

% D=sym('D',[3,3],'real');
% d=sym('d',[3,3],'real');
% lambda  =sym('lambda','real');
% g=sym('g','positive');
% 
% expression= det((D-lambda*eye(3,3)^2-(1/g_2)*d*d'));
% collected= collect(expression,lambda);
% lambda_solutions = solve(collected == 0, lambda);


D11=D(1,1);
D12=D(1,2);
D13=D(1,3);

D21=D(2,1);
D22=D(2,2);
D23=D(2,3);

D31=D(3,1);
D32=D(3,2);
D33=D(3,3);

d1=d(1,1);
d2=d(2,1);
d3=d(3,1);

D11_2=D11^2;
D12_2=D12^2;
D13_2=D13^2;

D21_2=D21^2;
D22_2=D22^2;
D23_2=D23^2;

D31_2=D31^2;
D32_2=D32^2;
D33_2=D33^2;

d1_2=d1^2;
d2_2=d2^2;
d3_2=d3^2;

g_2=g^2;


a6=1;
a5= - ((2*D11*g_2 + 2*D22*g_2 + 2*D33*g_2))/g_2;
a4=- ((d1_2 + d2_2 + d3_2 - D11_2*g_2 - D22_2*g_2 - D33_2*g_2 - 4*D11*D22*g_2 + 2*D12*D21*g_2 - 4*D11*D33*g_2 + 2*D13*D31*g_2 - 4*D22*D33*g_2 + 2*D23*D32*g_2))/g_2 ; 
a3= + ((2*D11*d2_2 + 2*D11*d3_2 + 2*D22*d1_2 + 2*D22*d3_2 + 2*D33*d1_2 + 2*D33*d2_2 - 2*D11*D22_2*g_2 - 2*D11_2*D22*g_2 - 2*D11*D33_2*g_2 - 2*D11_2*D33*g_2 - 2*D22*D33_2*g_2 - 2*D22_2*D33*g_2 - 2*D12*d1*d2 - 2*D13*d1*d3 - 2*D21*d1*d2 - 2*D23*d2*d3 - 2*D31*d1*d3 - 2*D32*d2*d3 + 2*D11*D12*D21*g_2 + 2*D11*D13*D31*g_2 + 2*D12*D21*D22*g_2 - 8*D11*D22*D33*g_2 + 4*D11*D23*D32*g_2 + 4*D12*D21*D33*g_2 - 2*D12*D23*D31*g_2 - 2*D13*D21*D32*g_2 + 4*D13*D22*D31*g_2 + 2*D13*D31*D33*g_2 + 2*D22*D23*D32*g_2 + 2*D23*D32*D33*g_2))/g_2 ;
a2=+ ((- D11_2*d2_2 - D11_2*d3_2 - D22_2*d1_2 - D22_2*d3_2 - D33_2*d1_2 - D33_2*d2_2 + D11_2*D22_2*g_2 + D12_2*D21_2*g_2 + D11_2*D33_2*g_2 + D13_2*D31_2*g_2 + D22_2*D33_2*g_2 + D23_2*D32_2*g_2 - D12*D21*d1_2 - D12*D21*d2_2 - 4*D11*D22*d3_2 + 2*D12*D21*d3_2 - D13*D31*d1_2 - 4*D11*D33*d2_2 + 2*D13*D31*d2_2 - D13*D31*d3_2 - 4*D22*D33*d1_2 + 2*D23*D32*d1_2 - D23*D32*d2_2 - D23*D32*d3_2 + 4*D11*D22*D33_2*g_2 + 4*D11*D22_2*D33*g_2 - 2*D12*D21*D33_2*g_2 - 2*D13*D22_2*D31*g_2 + 4*D11_2*D22*D33*g_2 - 2*D11_2*D23*D32*g_2 + D11*D12*d1*d2 + D11*D13*d1*d3 + D11*D21*d1*d2 + D12*D22*d1*d2 + 4*D11*D23*d2*d3 - 3*D12*D23*d1*d3 - 3*D13*D21*d2*d3 + 4*D13*D22*d1*d3 + D11*D31*d1*d3 + D21*D22*d1*d2 + 4*D11*D32*d2*d3 - 3*D12*D31*d2*d3 + 4*D12*D33*d1*d2 - 3*D13*D32*d1*d2 + D13*D33*d1*d3 + D22*D23*d2*d3 - 3*D21*D32*d1*d3 + 4*D21*D33*d1*d2 + 4*D22*D31*d1*d3 - 3*D23*D31*d1*d2 + D22*D32*d2*d3 + D23*D33*d2*d3 + D31*D33*d1*d3 + D32*D33*d2*d3 - 2*D11*D12*D21*D22*g_2 - 4*D11*D12*D21*D33*g_2 + 2*D11*D12*D23*D31*g_2 + 2*D11*D13*D21*D32*g_2 - 4*D11*D13*D22*D31*g_2 + 2*D12*D13*D21*D31*g_2 - 2*D11*D13*D31*D33*g_2 - 4*D11*D22*D23*D32*g_2 - 4*D12*D21*D22*D33*g_2 + 2*D12*D21*D23*D32*g_2 + 2*D12*D22*D23*D31*g_2 + 2*D13*D21*D22*D32*g_2 - 4*D11*D23*D32*D33*g_2 + 2*D12*D23*D31*D33*g_2 + 2*D13*D21*D32*D33*g_2 - 4*D13*D22*D31*D33*g_2 + 2*D13*D23*D31*D32*g_2 - 2*D22*D23*D32*D33*g_2))/g_2 ;
a1=- ((- 2*D11*D22_2*d3_2 - 2*D11_2*D22*d3_2 - 2*D11*D33_2*d2_2 - 2*D11_2*D33*d2_2 - 2*D22*D33_2*d1_2 - 2*D22_2*D33*d1_2 + 2*D13*D22_2*d1*d3 + 2*D11_2*D23*d2*d3 + 2*D12*D33_2*d1*d2 + 2*D11_2*D32*d2*d3 + 2*D21*D33_2*d1*d2 + 2*D22_2*D31*d1*d3 + 2*D11*D22_2*D33_2*g_2 + 2*D11*D23_2*D32_2*g_2 + 2*D11_2*D22*D33_2*g_2 + 2*D11_2*D22_2*D33*g_2 + 2*D12_2*D21_2*D33*g_2 + 2*D13_2*D22*D31_2*g_2 + 2*D11*D12*D21*d3_2 + 2*D11*D13*D31*d2_2 + 2*D12*D21*D22*d3_2 - 2*D12*D21*D33*d1_2 + 2*D12*D23*D31*d1_2 + 2*D13*D21*D32*d1_2 - 2*D13*D22*D31*d1_2 - 2*D11*D23*D32*d2_2 - 2*D12*D21*D33*d2_2 + 2*D12*D23*D31*d2_2 + 2*D13*D21*D32*d2_2 - 2*D11*D23*D32*d3_2 + 2*D12*D23*D31*d3_2 + 2*D13*D21*D32*d3_2 - 2*D13*D22*D31*d3_2 + 2*D22*D23*D32*d1_2 + 2*D13*D31*D33*d2_2 + 2*D23*D32*D33*d1_2 - 2*D11*D12*D21*D33_2*g_2 - 2*D11*D13*D22_2*D31*g_2 - 2*D12*D13*D21_2*D32*g_2 - 2*D12*D13*D23*D31_2*g_2 - 2*D12_2*D21*D23*D31*g_2 - 2*D12*D21*D22*D33_2*g_2 - 2*D11_2*D22*D23*D32*g_2 - 2*D13*D21*D23*D32_2*g_2 - 2*D13_2*D21*D31*D32*g_2 - 2*D12*D23_2*D31*D32*g_2 - 2*D13*D22_2*D31*D33*g_2 - 2*D11_2*D23*D32*D33*g_2 - 2*D11*D12*D23*d1*d3 - 2*D11*D13*D21*d2*d3 + 2*D11*D13*D22*d1*d3 - 2*D11*D12*D31*d2*d3 + 2*D11*D12*D33*d1*d2 - 2*D11*D13*D32*d1*d2 + 2*D11*D22*D23*d2*d3 - 2*D12*D22*D23*d1*d3 - 2*D13*D21*D22*d2*d3 - 2*D11*D21*D32*d1*d3 + 2*D11*D21*D33*d1*d2 + 2*D11*D22*D31*d1*d3 - 2*D11*D23*D31*d1*d2 + 2*D11*D22*D32*d2*d3 - 2*D12*D22*D31*d2*d3 + 2*D12*D22*D33*d1*d2 - 2*D13*D22*D32*d1*d2 + 2*D11*D23*D33*d2*d3 - 2*D12*D23*D33*d1*d3 - 2*D13*D21*D33*d2*d3 + 2*D13*D22*D33*d1*d3 - 2*D21*D22*D32*d1*d3 + 2*D21*D22*D33*d1*d2 - 2*D22*D23*D31*d1*d2 + 2*D11*D32*D33*d2*d3 - 2*D12*D31*D33*d2*d3 - 2*D13*D32*D33*d1*d2 - 2*D21*D32*D33*d1*d3 + 2*D22*D31*D33*d1*d3 - 2*D23*D31*D33*d1*d2 - 4*D11*D12*D21*D22*D33*g_2 + 2*D11*D12*D21*D23*D32*g_2 + 2*D11*D12*D22*D23*D31*g_2 + 2*D11*D13*D21*D22*D32*g_2 + 2*D12*D13*D21*D22*D31*g_2 + 2*D11*D12*D23*D31*D33*g_2 + 2*D11*D13*D21*D32*D33*g_2 - 4*D11*D13*D22*D31*D33*g_2 + 2*D11*D13*D23*D31*D32*g_2 + 2*D12*D13*D21*D31*D33*g_2 - 4*D11*D22*D23*D32*D33*g_2 + 2*D12*D21*D23*D32*D33*g_2 + 2*D12*D22*D23*D31*D33*g_2 + 2*D13*D21*D22*D32*D33*g_2 + 2*D13*D22*D23*D31*D32*g_2))/g_2 ;
a0= - (D11_2*D22_2*d3_2 + D12_2*D21_2*d3_2 + D11_2*D33_2*d2_2 + D13_2*D31_2*d2_2 + D22_2*D33_2*d1_2 + D23_2*D32_2*d1_2 - D11_2*D22_2*D33_2*g_2 - D11_2*D23_2*D32_2*g_2 - D12_2*D21_2*D33_2*g_2 - D12_2*D23_2*D31_2*g_2 - D13_2*D21_2*D32_2*g_2 - D13_2*D22_2*D31_2*g_2 + D12*D21*D33_2*d1_2 + D13*D22_2*D31*d1_2 + D12*D21*D33_2*d2_2 + D11_2*D23*D32*d2_2 + D13*D22_2*D31*d3_2 + D11_2*D23*D32*d3_2 - 2*D11*D12*D21*D22*d3_2 - D11*D12*D23*D31*d2_2 - D11*D13*D21*D32*d2_2 + D12*D13*D21*D31*d2_2 - D11*D12*D23*D31*d3_2 - D11*D13*D21*D32*d3_2 + D12*D13*D21*D31*d3_2 + D12*D21*D23*D32*d1_2 - D12*D22*D23*D31*d1_2 - D13*D21*D22*D32*d1_2 - 2*D11*D13*D31*D33*d2_2 + D12*D21*D23*D32*d3_2 - D12*D22*D23*D31*d3_2 - D13*D21*D22*D32*d3_2 - D12*D23*D31*D33*d1_2 - D13*D21*D32*D33*d1_2 + D13*D23*D31*D32*d1_2 - D12*D23*D31*D33*d2_2 - D13*D21*D32*D33*d2_2 + D13*D23*D31*D32*d2_2 - 2*D22*D23*D32*D33*d1_2 - D11*D13*D22_2*d1*d3 - D12*D13*D21_2*d2*d3 - D11*D12*D33_2*d1*d2 - D12_2*D21*D23*d1*d3 - D12*D13*D31_2*d2*d3 - D11_2*D22*D23*d2*d3 - D11*D21*D33_2*d1*d2 - D11*D22_2*D31*d1*d3 - D12*D21_2*D32*d1*d3 + D13*D21_2*D32*d1*d2 - D12_2*D21*D31*d2*d3 + D12_2*D23*D31*d1*d2 - D12*D22*D33_2*d1*d2 + D12*D23*D31_2*d1*d3 - D13*D23*D31_2*d1*d2 - D11_2*D22*D32*d2*d3 - D13_2*D21*D31*d2*d3 + D13_2*D21*D32*d1*d3 + D12*D23_2*D31*d2*d3 - D12*D23_2*D32*d1*d3 + D13*D21*D32_2*d2*d3 - D13*D23*D32_2*d1*d2 - D13*D22_2*D33*d1*d3 - D11_2*D23*D33*d2*d3 - D21*D22*D33_2*d1*d2 - D13_2*D31*D32*d1*d2 - D21*D23*D32_2*d1*d3 - D11_2*D32*D33*d2*d3 - D23_2*D31*D32*d1*d2 - D22_2*D31*D33*d1*d3 + D11*D12*D21*D23*d2*d3 + D11*D12*D22*D23*d1*d3 + D11*D13*D21*D22*d2*d3 + D12*D13*D21*D22*d1*d3 + D11*D12*D21*D32*d2*d3 + D11*D12*D22*D31*d2*d3 - D11*D12*D23*D32*d1*d2 + D11*D13*D22*D32*d1*d2 - D12*D13*D22*D31*d1*d2 + D11*D12*D23*D33*d1*d3 + D11*D13*D21*D33*d2*d3 + D11*D13*D23*D31*d2*d3 - D11*D13*D23*D32*d1*d3 - D12*D13*D21*D33*d1*d3 + D11*D21*D22*D32*d1*d3 - D11*D21*D23*D32*d1*d2 + D11*D22*D23*D31*d1*d2 + D12*D21*D22*D31*d1*d3 - D13*D21*D22*D31*d1*d2 + D11*D12*D31*D33*d2*d3 + D11*D13*D31*D32*d2*d3 + D11*D13*D32*D33*d1*d2 + D12*D13*D31*D33*d1*d2 - D12*D21*D23*D33*d2*d3 + D12*D22*D23*D33*d1*d3 + D13*D21*D22*D33*d2*d3 - D13*D22*D23*D31*d2*d3 + D13*D22*D23*D32*d1*d3 + D11*D21*D32*D33*d1*d3 - D11*D23*D31*D32*d1*d3 + D11*D23*D31*D33*d1*d2 - D12*D21*D31*D33*d1*d3 + D13*D21*D31*D33*d1*d2 - D12*D21*D32*D33*d2*d3 + D12*D22*D31*D33*d2*d3 + D12*D23*D32*D33*d1*d2 - D13*D22*D31*D32*d2*d3 + D13*D22*D32*D33*d1*d2 + D21*D22*D32*D33*d1*d3 + D21*D23*D32*D33*d1*d2 + D22*D23*D31*D32*d1*d3 + D22*D23*D31*D33*d1*d2 + 2*D11*D12*D21*D22*D33_2*g_2 + 2*D11*D13*D21*D23*D32_2*g_2 + 2*D12*D13*D22*D23*D31_2*g_2 + 2*D11*D12*D23_2*D31*D32*g_2 + 2*D11*D13*D22_2*D31*D33*g_2 + 2*D12*D13*D21_2*D32*D33*g_2 + 2*D13_2*D21*D22*D31*D32*g_2 + 2*D12_2*D21*D23*D31*D33*g_2 + 2*D11_2*D22*D23*D32*D33*g_2 - 2*D11*D12*D21*D23*D32*D33*g_2 - 2*D11*D12*D22*D23*D31*D33*g_2 - 2*D11*D13*D21*D22*D32*D33*g_2 - 2*D11*D13*D22*D23*D31*D32*g_2 - 2*D12*D13*D21*D22*D31*D33*g_2 - 2*D12*D13*D21*D23*D31*D32*g_2)/g_2 ;


coeffs = [a6, a5, a4, a3, a2, a1, a0];

roots_result = roots(coeffs);


real_result=[];

for i=1:length(roots_result)

    if imag(roots_result(i))==0

        real_result=[real_result;real(roots_result(i))];

    end

end

result=min(real_result);

end



 







