function fval=functionROSObject(x)
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
fval=sum((100*(x(2:4)-x(1:3).^2).^2-(x(1:3)-1).^2).^2);
end