function fval=functionPVDObject(x)
% Pressure vessel design (PVD) problem
% variable_number=4;
% object_function=@(x) functionPVDObject(x);
% object_function_low=@(x) functionPVDObjectLow(x);
% A=[-1,0,0.0193,0;
%     0,-1,0.00954,0;];
% B=[0;0];
% Aeq=[];
% Beq=[];
% low_bou=[1.0;0.625;25;25];
% up_bou=[1.375;1;150;240];
% nonlcon_function=@(x) functionPVDNonlcon(x);
% cheapcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq);
% model_function=[];
%
% x_best=[1.3006;0.6429;67.3860;25.0000]; fval_best=8949.4;
%
x1=x(1);x2=x(2);x3=x(3);x4=x(4);
fval=0.6224*x1*x3*x4+1.7781*x2*x3^2+3.1661*x1^2*x4+19.84*x1^2*x3;
end