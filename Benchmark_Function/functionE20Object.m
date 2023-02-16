function fval=functionE20Object(x)
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
ss=[1.8;0.4;2;1.2;1.4;0.6;1.6;0.2;0.8;1;1.3;1.1;2;1.4;0.5;0.3;1.6;0.7;0.3;1.9];
sh=[0.3;0.4;0.2;0.6;1;0.9;0.2;0.8;0.5;0.7;0.4;0.3;0.7;1;0.9;0.6;0.2;0.8;0.2;0.5];
fval=sum(linspace(1,20,20)'.*sh.*(x-ss).^2);
end