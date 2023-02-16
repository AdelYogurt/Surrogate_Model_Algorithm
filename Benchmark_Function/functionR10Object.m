function fval=functionR10Object(x)
% 10D Rosenbrock function
% variable_number=10;
% object_function=@(x) functionR10Object(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=zeros(variable_number,1);
% up_bou=ones(variable_number,1);
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];
%
% x_best=[1,1,1,1,1,1,1,1,1,1] fval_best=0
%
fval=sum(100*(x(2:10)-x(1:9).^2).^2+(x(1:9)-1).^2);
end