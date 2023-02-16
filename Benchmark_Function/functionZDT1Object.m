function f = functionZDT1Object(x)
% ZDT problem
% variable_number is 30
%
% object_function=@(x) functionZDT1Object(x);
% variable_number=10;
% low_bou=zeros(10,1);
% up_bou=ones(10,1);
% nonlcon_function=[];
% cheapcon_function=[];
%
% pareto_front: f2=1-f1.^0.5
%
variable_number=10;
f=zeros(2,1);
f(1) = x(1);
g=1+9*(sum(x(2:variable_number))/(variable_number-1));
f(2) = g * (1 - sqrt(x(1) / g));
end