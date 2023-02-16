function fval=functionA20Object(x)
% variable_number=20;
% object_function_H=@(x) functionA20Object(x);
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
s=[1.3;0.1;1.4;0.8;1.7;1;1.5;0.6;2;0.4;1.3;0.3;1.5;0.9;1.9;1.1;1.7;0.7;2.1;0.5];
fval=-20*exp(-0.2*sqrt(sum((x-s).^2)/10))-exp(sum(cos(2*1.3*pi*(x-s))/10))+20+exp(1);
end