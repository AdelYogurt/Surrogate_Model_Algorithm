function fval=functionA10Object(x)
% variable_number=10;
% object_function=@(x) functionA10Object(x);
% object_function_low=@(x) functionA10ObjectLow(x);
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
x=x(:);
s=[1.3;0.1;1.4;0.8;1.7;1;1.5;0.6;2;0.4];
fval=-20*exp(-0.2*sqrt(sum((x-s).^2)/10))-exp(sum(cos(2*1.3*pi*(x-s))/10))+20+exp(1);
end