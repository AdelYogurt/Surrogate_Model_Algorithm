function fval=functionA20ObjectLow(x)
% variable_number=20
% object_function=@(x) functionA20Object(x);
% object_function_low=@(x) functionA20ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=zeros(variable_number,1);
% up_bou=ones(variable_number,1);
% nonlcon_function=[];
% cheapcon_function=[];
%
% fval_min=0
% 
fval=-20*exp(-0.2*sqrt(sum(x.^2)/10))-exp(sum(cos(2*pi*x)/10))+20+exp(1);
end