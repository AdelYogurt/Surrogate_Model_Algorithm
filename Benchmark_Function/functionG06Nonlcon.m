function [con,coneq]=functionG06Nonlcon(x)
% variable_number=2;
% object_function=@functionG06Object;
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[13;0];
% up_bou=[100;100];
% nonlcon_function=@functionG06Nonlcon;
% cheapcon_function=[];
% model_function=[];
% 
% x_best=[14.0950;0.8430] fval_best=-6.9618e+03
%
g1=-(x(1)-5)^2-(x(2)-5)^2+100;
g2=(x(1)-6)^2+(x(2)-5)^2-82.81;
con=[g1;g2];
coneq=[];
end