function fval=functionDP20ObjectLow(x)
% variable_number=20;
% object_function=@(x) functionDP20Object(x);
% object_function_low=@(x) functionDP20ObjectLow(x);
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
fval=(x(1)-1)^2+sum(linspace(2,20,19)'.*(2*x(2:end).^2-x(1:end-1)).^2);
end