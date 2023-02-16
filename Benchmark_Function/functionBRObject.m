function fval=functionBRObject(x)
% Branin function
% variable_number=2;
% object_function=@(x) functionBRObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-5;10];
% up_bou=[0;15];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];
%
% x_best=[-3.1416;12.2750] fval_best=0.3979
%
x1=x(1);
x2=x(2);
fval=(x2-5.1/4/pi/pi*x1^2+5/pi*x1-6)^2+10*(1-1/8/pi)*cos(x1)+10;
end