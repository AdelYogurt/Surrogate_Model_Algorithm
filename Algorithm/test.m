clc;
clear;
close all hidden;

benchmark = BenchmarkFunction();

variable_number = 30;
object_function = @(x) benchmark.singleAckley30Object(x);
A = [];
B = [];
Aeq = [];
Beq = [];
low_bou = -15*ones(1,variable_number);
up_bou = 20*ones(1,variable_number);
nonlcon_function = [];
cheapcon_function = [];

[x_best,fval_best,exitflag,output] = ga(object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,optimoptions('ga','HybridFcn','fmincon'))