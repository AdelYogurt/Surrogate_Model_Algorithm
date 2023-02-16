function fval=functionE20ObjectLow(x)
% variable_number=20;
% object_function=@(x) functionE20Object(x);
% object_function_low=@(x) functionE20ObjectLow(x);
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
fval=sum(linspace(1,20,20)'.*x.^2);
end