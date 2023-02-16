function fval=functionCOSObject(x)
x=x(:);
par=[2,-2];
fval=sum(sum((x-par').^2-10*cos(2*pi*x)+10));
end