function fval=functionPKObject(x)
%
% variable_number=2;
% object_function=@(x) functionPKObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-3;-3];
% up_bou=[3;3];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];
%
% x_best=[0.2283;-1.6255] fval_best=-6.5511
%
x1=x(1);
x2=x(2);
fval=3*(1-x1).^2.*exp(-(x1.^2) - (x2+1).^2) ... 
   - 10*(x1/5 - x1.^3 - x2.^5).*exp(-x1.^2-x2.^2) ... 
   - 1/3*exp(-(x1+1).^2 - x2.^2) ;
end