function fval=functionDP20Object(x)
% variable_number=20;
% object_function=@(x) functionDP20Object(x);
% object_function_low=@(x) functionDP20ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=ones(variable_number,1)*-30;
% up_bou=ones(variable_number,1)*30;
% nonlcon_function=[];
% cheapcon_function=[];
%
% fval_min=0
% 
s=[1.8;0.5;2;1.2;0.4;0.2;1.4;0.3;1.6;0.6;0.8;1;1.3;1.9;0.7;1.6;0.3;1.1;2;1.4];
fval=func_G2_DP20_L(x-s);
end