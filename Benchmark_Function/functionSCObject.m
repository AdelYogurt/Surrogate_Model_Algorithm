function fval=functionSCObject(x)
% Six-hump came back function
% variable_number=2;
% object_function=@(x) functionSCObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-2;-2];
% up_bou=[2;2];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];
%
% x_best=[0.0898;-0.7127] fval_best=-1.0316
%
x1=x(1);
x2=x(2);
fval=4*x1^2-2.1*x1^4+x1^6/3+x1*x2-4*x2^2+4*x2^4;
end