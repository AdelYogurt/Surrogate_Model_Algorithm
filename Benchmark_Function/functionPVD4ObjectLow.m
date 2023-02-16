function fval=functionPVD4ObjectLow(x)
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
fval=functionPVDObject(x)*0.9+0.5;
end