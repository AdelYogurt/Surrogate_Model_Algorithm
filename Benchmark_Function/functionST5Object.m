function fval=functionST5Object(x)
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
s=[0.28;0.59;0.47;0.16;0.32];
fval=func_H2_ST5_L(x-s);
end