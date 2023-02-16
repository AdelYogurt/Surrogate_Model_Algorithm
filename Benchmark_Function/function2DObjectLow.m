function fval=function2DObjectLow(x)
% Rosenbrock problem, x_number=2
% low_bou=[-5;-5], up_bou=[5;5]
% [-3.6483;-0.0685],272.5563
% 
x1=x(1);
x2=x(2);
fval=function2DObject([0.9*x1;0.8*x2])-(x1+1)^2;
end