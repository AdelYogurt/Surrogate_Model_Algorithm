function fval=functionRSObject(x)
% Rastrigin function
% variable_number=2;
% object_function=@(x) functionRSObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-1;-1];
% up_bou=[1;1];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];
%
% x_best=[0;0] fval_best=-2
%
x1=x(1);
x2=x(2);
fval=x1^2+x2^2-cos(18*x1)-cos(18*x2);
end