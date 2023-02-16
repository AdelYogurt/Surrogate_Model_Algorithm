function fval=function2DObject(x)
% Rosenbrock problem, x_number=2
% low_bou=[-5;-5], up_bou=[5;5]
% [-3.6483;-0.0685],272.5563
% 
x1=x(1);
x2=x(2);
fval=(x1^2+x2-11)^2+(x2^2+x1+20)^2;
end