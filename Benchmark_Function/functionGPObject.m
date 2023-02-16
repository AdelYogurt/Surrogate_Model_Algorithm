function fval=functionGPObject(x)
% variable_number=2;
% object_function=@functionGPObject;
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-2;-2];
% up_bou=[2;2];
% nonlcon_function=[];
% cheapcon_function=[];
%
% x_best=[0;-1] fval_best=3
%
x=x(:);
fval=0;
x1=x(1);
x2=x(2);
fval=(1+(x1+x2+1)^2*...
    (19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2))*...
    (30+(2*x1-3*x2)^2*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2));
end