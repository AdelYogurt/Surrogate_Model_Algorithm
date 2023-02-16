function [con,coneq]=functionPVD4Nonlcon(x)
% Pressure vessel design (PVD4) problem
% variable_number=4;
% object_function=@(x) functionPVD4Object(x);
% object_function_low=@(x) functionPVD4ObjectLow(x);
% A=[-1,0,0.0193,0;
%     0,-1,0.00954,0;];
% B=[0;0];
% Aeq=[];
% Beq=[];
% low_bou=[0;0;0;0];
% up_bou=[1;1;50;240];
% nonlcon_function=@(x) functionPVD4Nonlcon(x);
% cheapcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq);
% model_function=[];
%
% x_best=[0.7276;0.3596;37.6991;240.0000]; fval_best=5804.45;
%
x1=x(1);x2=x(2);x3=x(3);x4=x(4);
g3=-pi*x3^2*x4-4/3*pi*x3^3+1296000;
if g3 >=0
    g3=log(1+g3);
else
    g3=-log(1-g3);
end
con=g3;
coneq=[];
end