function fval=functionEP20ObjectLow(x)
% variable_number=20;
% object_function=@(x) functionEP20Object(x);
% object_function_low=@(x) functionEP20ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=ones(1,variable_number)*-30;
% up_bou=ones(1,variable_number)*30;
% nonlcon_function=[];
% cheapcon_function=[];
%
% fval_min=0
%
fval=sum(linspace(1,20,20)'.*x.^2);
end