function fval=functionCOLObjectLow(x)
% variable_number=4;
% object_function=@(x) functionCOLObject(x);
% object_function_low=@(x) functionCOLObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[0;0;0;0];
% up_bou=[1;1;1;1];
% nonlcon_function=[];
% cheapcon_function=[];
% 
fval=func_H2_COL_H([0.8;0.8;0.5;0.5].*x);
end