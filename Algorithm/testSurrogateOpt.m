clc;
clear;
close all hidden;

benchmark = BenchmarkFunction();

% variable_number = 2;
% object_function = @(x) benchmark.singleGPObject(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = [-2, -2];
% up_bou = [2, 2];
% nonlcon_function = [];
% cheapcon_function = [];

% variable_number = 20;
% object_function = @(x) benchmark.singleEP20Object(x);
% object_function_LF = @(x) benchmark.singleEP20ObjectLow(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = ones(1, variable_number)*-30;
% up_bou = ones(1, variable_number)*30;
% nonlcon_function = [];
% nonlcon_function_LF = [];
% cheapcon_function = [];

% variable_number = 6;
% object_function = @(x) benchmark.singleHNObject(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = zeros(1,variable_number);
% up_bou = ones(1,variable_number);
% nonlcon_function = [];
% cheapcon_function = [];
% model_function = [];

% variable_number = 30;
% object_function = @(x) sum((x-rand(1,30)).^2,2);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = ones(1, variable_number)*-30;
% up_bou = ones(1, variable_number)*30;
% nonlcon_function = [];
% cheapcon_function = [];

variable_number = 13;
object_function = @(x) benchmark.singleG01Object(x);
object_function_low = @(x) benchmark.singleG01ObjectLow(x);
A = [ 2   2   0   0   0   0   0   0   0   1   1   0   0;
    2   0   2   0   0   0   0   0   0   1   0   1   0;
    0   2   2   0   0   0   0   0   0   0   1   1   0;
    -8  0   0   0   0   0   0   0   0   1   0   0   0;
    0   -8  0   0   0   0   0   0   0   0   1   0   0;
    0   0   0   -2  -1  0   0   0   0   1   0   0   0;
    0   0   0   0   0   -2  -1  0   0   0   1   0   0;
    0   0   0   0   0   0   0   -2  -1  0   0   1   0;
    ];
B = [10;10;10;0;0;0;0;0];
Aeq = [];
Beq = [];
low_bou = zeros(1,13);
up_bou = ones(1,13);
up_bou(10:12) = 100;
nonlcon_function = @(x) cheapconFunction(x,A,B,Aeq,Beq,[]);
nonlcon_function_LF = @(x) cheapconFunction(x,A,B,Aeq,Beq,[]);
cheapcon_function = [];

repeat_number = 10;
result_fval = zeros(repeat_number, 1);
max_NFE = 200;

objconstr = @(x) modelFunction(x,object_function,nonlcon_function);

for repeat_index = 1:repeat_number
    x_initial_list = lhsdesign(2*variable_number,variable_number,'iterations',100,'criterion','maximin').*(up_bou-low_bou)+low_bou;
    surrogateopt_option = optimoptions('surrogateopt', 'MaxFunctionEvaluations', max_NFE, 'Display', 'none','PlotFcn','surrogateoptplot','InitialPoints',x_initial_list);

    [x_best, fval_best, exitflag, output] = surrogateopt...
        (objconstr, low_bou, up_bou, surrogateopt_option);

    result_fval(repeat_index) = fval_best;
end

fprintf('Fval     : lowest = %4.4f, mean = %4.4f, worst = %4.4f, std = %4.4f \n', min(result_fval), mean(result_fval), max(result_fval), std(result_fval));
object_function_name=char(object_function);
save([object_function_name(15:end-3),'_',num2str(max_NFE),'_surrogatopt','.mat']);


function objcon = modelFunction(x,object_function,nonlcon_function)
fval = object_function(x);
if isempty(nonlcon_function)
    objcon = fval;
else
    [con,coneq] = nonlcon_function(x);
    objcon.Fval = fval;
    objcon.Ineq = con;
end
end
