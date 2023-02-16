function fval=functionTNKObject(x)
% TNK problem
% variable_number is 2
%
% object_function=@(x) functionTNKObject(x);
% variable_number=2;
% low_bou=zeros(1,2);
% up_bou=ones(1,2)*pi;
% nonlcon_function=@(x) functionTNKNonlcon(x);
% cheapcon_function=[];
%
fval=zeros(2,1);
fval(1)=x(1);
fval(2)=x(2);
end