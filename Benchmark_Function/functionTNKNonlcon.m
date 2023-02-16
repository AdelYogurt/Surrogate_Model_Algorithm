function [con,coneq]=functionTNKNonlcon(x)
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
con=zeros(2,1);
x1=x(1);
x2=x(2);
if x1 == 0 && x2 == 0
    con(1)=-(x1^2+x2^2-1-0.1);
else
    con(1)=-(x1^2+x2^2-1-0.1*cos(16*atan(x1/x2)));
end
con(2)=(x1-0.5)^2+(x2-0.5)^2-0.5;
con=con*100;
coneq=[];
end