function fval=functionROSObjectLow(x)
% Rosenbrock problem
% variable_number=4
% object_function=@(x) functionROSObject(x);
% object_function_low=@(x) functionROSObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[2;2;2;2]
% up_bou=[2;2;2;2]
% nonlcon_function=[];
% cheapcon_function=[];
% 
% fval_min=0
% 
fval=sum((100*(0.5*x(2:4)-0.6*x(1:3).^2).^2-(0.5*x(1:3)-0.5).^2).^2);
end