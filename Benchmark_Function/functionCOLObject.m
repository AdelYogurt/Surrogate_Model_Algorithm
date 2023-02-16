function fval=functionCOLObject(x)
% variable_number=4;
% object_function=@(x) functionCOLObject(x);
% object_function_low=@(x) functionCOLObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[0;0;0;0];
% up_bou=[1;1;1;1];
% nonlcon_function=[];
% cheapcon_function=[];
%
% fval_min=-9.3361;
%

x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
fval=100*(x1^2-x2)^2+(x1-1)^2+(x3-1)^2+90*(x3^2-x4)^2+...
    10.1*((x2-1)^2-(x4-1)^2)+19.8*(x2-1)*(x4-1);
end