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

variable_number = 20;
object_function = @(x) benchmark.singleEP20Object(x);
object_function_LF = @(x) benchmark.singleEP20ObjectLow(x);
A = [];
B = [];
Aeq = [];
Beq = [];
low_bou = ones(1, variable_number)*-30;
up_bou = ones(1, variable_number)*30;
nonlcon_function = [];
nonlcon_function_LF = [];
cheapcon_function = [];

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

repeat_number = 10;
result_fval = zeros(repeat_number, 1);
max_NFE = 200;

for repeat_index = 1:repeat_number
    x_initial_list = lhsdesign(ceil(variable_number/10)*20,variable_number,'iterations',100,'criterion','maximin').*(up_bou-low_bou)+low_bou;
    surrogateopt_option = optimoptions('surrogateopt', 'MaxFunctionEvaluations', max_NFE, 'Display', 'none','PlotFcn','surrogateoptplot','InitialPoints',x_initial_list);

    [x_best, fval_best, exitflag, output] = surrogateopt...
        (object_function, low_bou, up_bou, surrogateopt_option);

    result_fval(repeat_index) = fval_best;
end

fprintf('Fval     : lowest = %4.4f, mean = %4.4f, worst = %4.4f, std = %4.4f \n', min(result_fval), mean(result_fval), max(result_fval), std(result_fval));
object_function_name=char(object_function);
save([object_function_name(15:end-3),'_',num2str(max_NFE),'_surrogatopt','.mat']);