function ab = quaternProd(a, b)

%  b=quaternProd([-a(1),a(2:4)], ab)

% ab(1)  =[b(1), -b(2), -b(3), -b(4);...   a(1)
% ab(2)    b(2),  b(1),  b(4), -b(3);...   a(2)
% ab(3)    b(3), -b(4),  b(1),  b(2);...   a(3)
% ab(4)    b(4),  b(3), -b(2),  b(1)]      a(4)

ab=a;


ab(1) = a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
ab(2) = a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
ab(3) = a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
ab(4) = a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1);
if ab(1)<0
    ab=-ab;
end


end







