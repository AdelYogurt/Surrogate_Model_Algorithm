function fval=functionGFObject(x)
x=x(:);
c1=1.5;
c2=2.25;
c3=2.625;
x1=x(1);
x2=x(2);
u1=c1-x1*(1-x2^1);
u2=c2-x1*(1-x2^2);
u3=c3-x1*(1-x2^3);
fval=u1^2+u2^2+u3^2;
end