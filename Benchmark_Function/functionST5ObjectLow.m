function fval=functionST5ObjectLow(x)
% variable_number=5;
% object_function=@(x) functionST5Object(x);
% object_function_low=@(x) functionST5ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-5;-5;-5;-5;-5];
% up_bou=[5;5;5;5;5];
% nonlcon_function=[];
% cheapcon_function=[];
%
% fval_min=0
% 
fval=0.5*sum(x.^4-16*x.^2+5*x);
end