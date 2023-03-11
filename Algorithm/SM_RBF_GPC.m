clc;
clear;
close all hidden;

data_library_name = 'optimal_data_library';

benchmark = BenchmarkFunction();

% variable_number = 2;
% object_function = @(x) benchmark.singleGPObject(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = [-2,-2];
% up_bou = [2,2];
% nonlcon_function = [];
% cheapcon_function = [];

% variable_number = 2;
% object_function = @(x) benchmark.singlePKObject(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = [-3,-3];
% up_bou = [3,3];
% nonlcon_function = [];
% cheapcon_function = [];

% variable_number = 2;
% object_function = @(x) benchmark.single2DObject(x);
% object_function_LF = @(x) benchmark.single2DObjectLow(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = [-5,-5];
% up_bou = [5,5];
% nonlcon_function = [];
% nonlcon_function_LF = [];
% cheapcon_function = [];

% variable_number = 4;
% object_function = @(x) benchmark.singleROSObject(x);
% object_function_LF = @(x) benchmark.singleROSObjectLow(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = [-2,-2,-2,-2];
% up_bou = [2,2,2,2];
% nonlcon_function = [];
% nonlcon_function_LF = [];
% cheapcon_function = [];

% variable_number = 20;
% object_function = @(x) benchmark.singleDP20Object(x);
% object_function_LF = @(x) benchmark.singleDP20ObjectLow(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = ones(1,variable_number)*-30;
% up_bou = ones(1,variable_number)*30;
% nonlcon_function = [];
% nonlcon_function_LF = [];
% cheapcon_function = [];

variable_number = 20;
object_function = @(x) benchmark.singleEP20Object(x);
object_function_LF = @(x) benchmark.singleEP20ObjectLow(x);
A = [];
B = [];
Aeq = [];
Beq = [];
low_bou = ones(1,variable_number)*-30;
up_bou = ones(1,variable_number)*30;
nonlcon_function = [];
nonlcon_function_LF = [];
cheapcon_function = [];

% variable_number = 30;
% object_function = @(x) benchmark.singleAckley30Object(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = -15*ones(1,variable_number);
% up_bou = 20*ones(1,variable_number);
% nonlcon_function = [];
% cheapcon_function = [];

% variable_number = 2;
% object_function = @(x) benchmark.singleG06Object(x);
% object_function_LF = @(x) benchmark.singleG06ObjectLow(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = [13,0];
% up_bou = [100,100];
% nonlcon_function = @(x) benchmark.singleG06Nonlcon(x);
% nonlcon_function_LF = @(x) benchmark.singleG06NonlconLow(x);
% cheapcon_function = [];
% model_function = [];

% variable_number = 4;
% object_function = @(x) benchmark.singlePVD4Object(x);
% A = [-1,0,0.0193,0;
%     0,-1,0.00954,0;];
% B = [0;0];
% Aeq = [];
% Beq = [];
% low_bou = [0,0,0,0];
% up_bou = [1,1,50,240];
% % nonlcon_function = @(x) benchmark.singlePVD4Nonlcon(x);
% nonlcon_function = @(x) cheapconFunction(x,A,B,Aeq,Beq,@(x) benchmark.singlePVD4Nonlcon(x));
% cheapcon_function = [];
% model_function = [];

% variable_number = 13;
% object_function = @(x) benchmark.singleG01Object(x);
% object_function_low = @(x) benchmark.singleG01ObjectLow(x);
% A = [
%     2   2   0   0   0   0   0   0   0   1   1   0   0;
%     2   0   2   0   0   0   0   0   0   1   0   1   0;
%     0   2   2   0   0   0   0   0   0   0   1   1   0;
%     -8  0   0   0   0   0   0   0   0   1   0   0   0;
%     0   -8  0   0   0   0   0   0   0   0   1   0   0;
%     0   0   -8  0   0   0   0   0   0   0   0   1   0
%     0   0   0   -2  -1  0   0   0   0   1   0   0   0;
%     0   0   0   0   0   -2  -1  0   0   0   1   0   0;
%     0   0   0   0   0   0   0   -2  -1  0   0   1   0;
%     ];
% B = [10;10;10;0;0;0;0;0;0];
% Aeq = [];
% Beq = [];
% low_bou = zeros(1,13);
% up_bou = ones(1,13);
% up_bou(10:12) = 100;
% nonlcon_function = @(x) cheapconFunction(x,A,B,Aeq,Beq,[]);
% nonlcon_function_LF = @(x) cheapconFunction(x,A,B,Aeq,Beq,[]);
% cheapcon_function = [];
 
% x_initial = rand(1,variable_number).*(up_bou-low_bou)+low_bou;
% [x_best,fval_best,~,output] = fmincon(object_function,x_initial,A,B,Aeq,Beq,low_bou,up_bou,[],optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',10000,'Display','iter-detailed'))

%% single run

delete([data_library_name,'.txt']);
delete('result_total.txt');

[x_best,fval_best,NFE,output] = optimalSurrogateRBFGPC...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,[],200,500)
result_x_best = output.result_x_best;
result_fval_best = output.result_fval_best;

figure(1);
plot(result_fval_best);

figure(2);
[x_list,fval_list,con_list,coneq_list] = dataLibraryLoad...
    (data_library_name,low_bou,up_bou);
scatter3(x_list(:,1),x_list(:,2),fval_list);
xlabel('X');
ylabel('Y');
zlabel('Z');

%% repeat run

% repeat_number = 10;
% result_fval = zeros(repeat_number,1);
% max_NFE = 200;
% for repeat_index = 1:repeat_number
%     delete([data_library_name,'.txt']);
%     delete('result_total.txt');
% 
%     [x_best,fval_best,NFE,output] = optimalSurrogateRBFGPC...
%         (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
%         cheapcon_function,[],max_NFE,300);
%     
%     result_fval(repeat_index) = fval_best;
% end
% 
% fprintf('Fval     : lowest = %4.4f,mean = %4.4f,worst = %4.4f,std = %4.4f \n',min(result_fval),mean(result_fval),max(result_fval),std(result_fval));
% object_function_name = char(object_function);
% save([object_function_name(15:end-3),'_',num2str(max_NFE),'_RBF_GPC','.mat']);

%% main
function [x_best,fval_best,NFE,output] = optimalSurrogateRBFGPC...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance,x_initial_list)
% KRG-CDE optimization algorithm
%
% referance: [1] 叶年辉,龙腾,武宇飞,et al.
% 基于Kriging代理模型的约束差分进化算法 [J]. 航空学报,2021,42(6): 13.
%
% Copyright 2022 Adel
%
if nargin < 12
    x_initial_list = [];
    if nargin < 11 || isempty(nonlcon_torlance)
        nonlcon_torlance = 1e-3;
        if nargin < 10 || isempty(torlance)
            torlance = 1e-3;
            if nargin < 9
                iteration_max = [];
                if nargin < 8
                    NFE_max = [];
                end
            end
        end
    end
end

if nargin < 7
    model_function = [];
    if nargin < 6
        cheapcon_function = [];
        if nargin < 5
            nonlcon_function = [];
        end
    end
end

DRAW_FIGURE_FLAG = 0; % whether draw data
INFORMATION_FLAG = 0; % whether print data
CONVERGENCE_JUDGMENT_FLAG = 0; % whether judgment convergence

% hyper parameter
sample_number_initial = 60;
trial_number = min(100*variable_number,100);
coord_select_prob_initial = min(20/variable_number,1);
sigma_coord_initial = 0.1*(up_bou-low_bou);
sigma_coord_min = 0.2*1/64*(up_bou-low_bou);
sigma_coord_max = 2*(up_bou-low_bou);
tau_success = 3;
tau_fail = max(variable_number,5);
RBF_number = max(100,(variable_number+1)*(variable_number+2)/2);
scaling_factor = 0.8; % F
cross_rate = 0.8;

% max fval when normalize fval,con,coneq
nomlz_fval = 10;

protect_range = 1e-6;

data_library_name = 'optimal_data_library';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done = 0;NFE = 0;iteration = 0;

% if do not input model_function,generate model_function
if isempty(model_function)
    model_function = @(x) modelFunction(x,object_function,nonlcon_function);
end

% step 2
% generate initial sample x_list
if isempty(x_initial_list)
%     [~,x_updata_list,~] = getLatinHypercube...
%         (sample_number_initial,variable_number,[],low_bou,up_bou,cheapcon_function);
    x_updata_list = lhsdesign(sample_number_initial,variable_number).*(up_bou-low_bou)+low_bou;
else
    x_updata_list = x_initial_list;
end

% detech expensive constraints
if ~isempty(x_updata_list)
    [~,con,coneq] = dataLibraryUpdata...
        (data_library_name,model_function,x_updata_list(1,:));NFE = NFE+1;
    x_updata_list = x_updata_list(2:end,:);
else
    [~,con,coneq] = dataLibraryLoad(data_library_name,low_bou,up_bou);
end
if ~isempty(con) || ~isempty(coneq)
    expensive_nonlcon_flag = 1;
else
    expensive_nonlcon_flag = 0;
end

% NFE and iteration setting
if isempty(NFE_max)
    if expensive_nonlcon_flag
        NFE_max = 50*variable_number;
    else
        NFE_max = 20*variable_number;
    end
end

if isempty(iteration_max)
    if expensive_nonlcon_flag
        iteration_max = 50*variable_number;
    else
        iteration_max = 20*variable_number;
    end
end
result_x_best = zeros(iteration_max,variable_number);
result_fval_best = zeros(iteration_max,1);

% import data from data library
[x_list,fval_list,con_list,coneq_list] = dataLibraryLoad...
    (data_library_name,low_bou,up_bou);

% updata data library by x_list
[fval_updata_list,con_updata_list,coneq_updata_list] = dataLibraryUpdata...
    (data_library_name,model_function,x_updata_list);NFE = NFE+size(x_updata_list,1);
x_list = [x_list;x_updata_list];
fval_list = [fval_list;fval_updata_list];
vio_list = zeros(size(x_list,1),1);
if ~isempty(con_list)
    con_list = [con_list;con_updata_list];
    vio_list = vio_list+sum(max(con_list-nonlcon_torlance,0),2);
end
if ~isempty(coneq_list)
    coneq_list = [coneq_list;coneq_updata_list];
    vio_list = vio_list+sum((abs(coneq_list)-nonlcon_torlance),2);
end

% step 0
% find best result to record
[x_best,fval_best,con_best,coneq_best] = findMinRaw...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance);

iteration = iteration+1;
sigma_coord = sigma_coord_initial;
C_success = 0;
C_fail = 0;
search_flag = 'G'; % 'G' is global search,'l' is local search
while ~done
    infor_search_flag = search_flag;
    % nomalization all data by max fval
    fval_max = max(abs(fval_list),[],1);
    fval_nomlz_list = fval_list./fval_max*nomlz_fval;
    if ~isempty(con_list)
        con_max_list = max(abs(con_list),[],1);
        con_nomlz_list = con_list./con_max_list*nomlz_fval;
    else
        con_nomlz_list = [];
    end
    if ~isempty(coneq_list)
        coneq_max_list = max(abs(coneq_list),[],1);
        coneq_nomlz_list = coneq_list./coneq_max_list*nomlz_fval;
    else
        coneq_nomlz_list = [];
    end

    % find fesiable data in current data library
    feasi_boolean_list = vio_list <= 0;

    % step 4
    % construct RBF model
    % base on distance to x_list to repace predict variance
    [RBF_model_fval,RBF_model_con,RBF_model_coneq,output_RBF] = getRadialBasisModel...
        (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list);
    object_function_surrogate = output_RBF.object_function_surrogate;
    nonlcon_function_surrogate = output_RBF.nonlcon_function_surrogate;

    improve_flag = 0;
%     search_flag = 'G';
    if search_flag == 'G'
        % global search

        % step 5
        % perturbed to get trial point
        coord_select_prob = coord_select_prob_initial*(1-log(NFE-sample_number_initial+1)/log(NFE_max-sample_number_initial));

        % select coordinates to perturb
        coordinate_select_list = [];
        for variable_index = 1:variable_number
            if rand() < coord_select_prob
                coordinate_select_list = [coordinate_select_list,variable_index];
            end
        end
        if isempty(coordinate_select_list)
            coordinate_select_list = randi(variable_number);
        end

        % generate trial point
        x_trial_list = repmat(x_best,trial_number,1);
        for coordinate_select_index = 1:length(coordinate_select_list)
            coordinate_select = coordinate_select_list(coordinate_select_index);
            x_trial_list(:,coordinate_select) = x_trial_list(:,coordinate_select)+...
                normrnd(0,sigma_coord(coordinate_select),[trial_number,1]);
        end
%         x_trial_list = repmat(x_best,trial_number,1);
%         for variable_index = 1:variable_number
%             x_trial_list(:,variable_index) = x_trial_list(:,variable_index)+...
%                 normrnd(0,sigma_coord(variable_index),[trial_number,1]);
%         end
        x_trial_list = max(x_trial_list,low_bou);
        x_trial_list = min(x_trial_list,up_bou);
       
        R2 = 1-sum((RBF_model_fval.beta./diag(RBF_model_fval.inv_radialbasis_matrix)).^2)/sum((RBF_model_fval.Y-mean(RBF_model_fval.Y)).^2);
        fprintf('global R2:  %f\n',R2);

        % evaluate each x_offspring fval and constraints
        [fval_pred_trial_list] = object_function_surrogate(x_trial_list);

        % variance were replaced by x distance to exist x
        % distance scale
        dis = zeros(size(x_trial_list,1),size(x_list,1));
        for vari_index = 1:variable_number
            dis = dis+(x_trial_list(:,vari_index)-x_list(:,vari_index)').^2;
        end
        D = min(sqrt(dis),[],2);
        D_min = min(D); D_max = max(D);
        fval_var_trial_list = D;
        if expensive_nonlcon_flag

            if ~isempty(nonlcon_function_surrogate)
                [con_pred_trial_list,coneq_pred_trial_list] = nonlcon_function_surrogate(x_trial_list);
                con_var_trial_list = fval_var_trial_list;
                coneq_var_trial_list = fval_var_trial_list;
            end

            vio_trial_list = zeros(trial_number,1);
            if ~isempty(con_nomlz_list)
                vio_trial_list = vio_trial_list+sum(max(con_pred_trial_list-nonlcon_torlance,0),2);
            end
            if ~isempty(coneq_nomlz_list)
                vio_trial_list = vio_trial_list+sum((abs(con_pred_trial_list)-nonlcon_torlance),2);
            end

            feasi_boolean_trial_list = vio_trial_list <= nonlcon_torlance;
        else
            feasi_boolean_trial_list = true(ones(1,trial_number));
        end
        
        % if have feasiable_index_list,only use feasiable to choose
        if all(~feasi_boolean_trial_list)
            % base on constaints improve select global infill
            % lack process of equal constraints
            con_nomlz_base = max(min(con_nomlz_list,[],1),0);
            con_impove_probability_list = sum(...
                normcdf((con_nomlz_base-con_pred_trial_list)./sqrt(con_var_trial_list)),2);
            [~,con_best_index] = max(con_impove_probability_list);
            con_best_index = con_best_index(1);
            x_global_infill = x_trial_list(con_best_index,:);
        else
            % base on fitness DE point to select global infill
            if expensive_nonlcon_flag
                x_trial_list = x_trial_list(feasi_boolean_trial_list,:);
                fval_pred_trial_list = fval_pred_trial_list(feasi_boolean_trial_list);
                fval_var_trial_list = fval_var_trial_list(feasi_boolean_trial_list);
            end

            fval_pred_trial_min = min(fval_pred_trial_list,[],1);
            fval_pred_trial_max = max(fval_pred_trial_list,[],1);
            fval_var_trial_min = min(fval_var_trial_list,[],1);
            fval_var_trial_max = max(fval_var_trial_list,[],1);
            % modify
            fitness_list = (fval_pred_trial_list-fval_pred_trial_min)/(fval_pred_trial_max-fval_pred_trial_min)+...
                (fval_var_trial_max-fval_var_trial_list)/(fval_var_trial_max-fval_var_trial_min);
            [~,fitness_best_index] = min(fitness_list);
            fitness_best_index = fitness_best_index(1);
            x_global_infill = x_trial_list(fitness_best_index,:);
        end
        
        [x_global_infill,fval_global_infill,con_global_infill,coneq_global_infill,NFE_p,repeat_index] = dataLibraryUpdataProtect...
            (data_library_name,model_function,x_global_infill,...
            x_list,low_bou,up_bou,protect_range);NFE = NFE+NFE_p;
        
        x_list = [x_list;x_global_infill];
        fval_list = [fval_list;fval_global_infill];
        vio_infill = 0;
        if ~isempty(con_list)
            con_list = [con_list;con_global_infill];
            vio_infill = vio_infill+sum(max(con_global_infill-nonlcon_torlance,0),2);
        end
        if ~isempty(coneq_list)
            coneq_list = [coneq_list;coneq_global_infill];
            vio_infill = vio_infill+sum((abs(coneq_global_infill)-nonlcon_torlance),2);
        end
        vio_list = [vio_list;vio_infill];

        % whether impove pupolation,if imporve,continue global
        % notice last one is x_local fval,con and coneq
        next_search_flag = 'l';
        if isempty(x_global_infill)
            % continue;
            x_global_infill = x_list(repeat_index,:);
            fval_global_infill = fval_list(repeat_index,:);
            if ~isempty(con_list)
                con_global_infill = con_list(repeat_index,:);
            end
            if ~isempty(coneq_list)
                coneq_global_infill = coneq_list(repeat_index,:);
            end
        end

        % infill point violation
        
        if all(~feasi_boolean_list)
            min_vio = min(vio_list);

            % improve, continue global search
            if vio_infill < min_vio
                next_search_flag = 'G';
                improve_flag = 1;
            end
        else
            if expensive_nonlcon_flag
                min_fval = min(fval_list([feasi_boolean_list;false(1)]),[],1);
                if vio_infill <= 0 % if exist constraint, infill point should feasible
                    vio_judge = 1;
                else
                    vio_judge = 0;
                end
            else
                min_fval = min(fval_list(1:end-1));
                vio_judge = 1;
            end

            % imporve, continue global search
            if fval_global_infill < min_fval && vio_judge
                next_search_flag = 'G';
                improve_flag = 1;
            end
        end
    
        if DRAW_FIGURE_FLAG && variable_number < 3
            interpVisualize(RBF_model_fval,low_bou,up_bou);
            line(x_global_infill(1),x_global_infill(2),fval_global_infill./fval_max*nomlz_fval,'Marker','o','color','r');
            hold on;
            scatter3(x_trial_list(:,1),x_trial_list(:,2),fitness_list*10);
            hold off;
        end
    else
        % local search
        [~,~,~,~,~,index_list] = rankData...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            cheapcon_function,nonlcon_torlance);

        % step 8
        % rand select initial local point from x_list
%         x_index = index_list(randi(min(size(x_list,1),2*variable_number)));
        x_index = index_list(1);
        x_initial = x_list(x_index,:);
        
        % get local infill point
        % obtian total constraint function
        if ~isempty(nonlcon_function_surrogate) || ~isempty(cheapcon_function)
            constraint_function = @(x) totalconFunction...
                (x,nonlcon_function_surrogate,cheapcon_function);
        else
            constraint_function = [];
        end
        fmincon_options = optimoptions('fmincon','Display','none','Algorithm','sqp','MaxIterations',50);
        x_local_infill = fmincon(object_function_surrogate,x_initial,[],[],[],[],...
            low_bou,up_bou,constraint_function,fmincon_options);
        
        [x_local_infill,fval_local_infill,con_local_infill,coneq_local_infill,NFE_p,repeat_index] = dataLibraryUpdataProtect...
            (data_library_name,model_function,x_local_infill,...
            x_list,low_bou,up_bou,protect_range);NFE = NFE+NFE_p;

        R2 = 1-sum((RBF_model_fval.beta./diag(RBF_model_fval.inv_radialbasis_matrix)).^2)/sum((RBF_model_fval.Y-mean(RBF_model_fval.Y)).^2);
        fprintf('local R2:  %f\n',R2);
        
        x_list = [x_list;x_local_infill];
        fval_list = [fval_list;fval_local_infill];
        vio_infill = 0;
        if ~isempty(con_list)
            con_list = [con_list;con_local_infill];
            vio_infill = vio_infill+sum(max(con_local_infill-nonlcon_torlance,0),2);
        end
        if ~isempty(coneq_list)
            coneq_list = [coneq_list;coneq_local_infill];
            vio_infill = vio_infill+sum((abs(coneq_local_infill)-nonlcon_torlance),2);
        end
        vio_list = [vio_list;vio_infill];
        
        % whether impove pupolation,if imporve,continue local
        % notice last one is x_local fval,con and coneq
        next_search_flag = 'G';
        if isempty(x_local_infill)
            % continue;
            x_local_infill = x_list(repeat_index,:);
            fval_local_infill = fval_list(repeat_index,:);
            if ~isempty(con_list)
                con_local_infill = con_list(repeat_index,:);
            end
            if ~isempty(coneq_list)
                coneq_local_infill = coneq_list(repeat_index,:);
            end
        end

        if all(~feasi_boolean_list)
            min_vio = min(vio_list);

            % improve, continue local search
            if vio_infill < min_vio
                next_search_flag = 'l';
                improve_flag = 1;
            end
        else
            if expensive_nonlcon_flag
                min_fval = min(fval_list([feasi_boolean_list;false(1)]),[],1);
                if vio_infill <= 0 % if exist constraint, infill point should feasible
                    vio_judge = 1;
                else
                    vio_judge = 0;
                end
            else
                min_fval = min(fval_list(1:end-1));
                vio_judge = 1;
            end

            % imporve, continue local search
            if fval_local_infill < min_fval && vio_judge
                next_search_flag = 'l';
                improve_flag = 1;
            end
        end

        if DRAW_FIGURE_FLAG && variable_number < 3
            interpVisualize(RBF_model_fval,low_bou,up_bou);
            line(x_local_infill(1),x_local_infill(2),fval_local_infill./fval_max*nomlz_fval,'Marker','o','color','r');
        end
    end
    
    % step 7
    % adjust step size
%     if ~improve_flag
%         C_success = C_success+1;
%         C_fail = 0;
%     else
%         C_success = 0;
%         C_fail = C_fail+1;
%     end
% 
%     if C_success >= tau_success
%         sigma_coord = min(2*sigma_coord,sigma_coord_max);
%         C_success = 0;
%     end
% 
%     if C_fail >= tau_fail
%         sigma_coord = max(sigma_coord/2,sigma_coord_min);
%         C_fail = 0;
%     end

    % find best result to record
    [x_best,fval_best,con_best,coneq_best] = findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);
    vio_best = 0;
    if ~isempty(con_list)
        vio_best = vio_best+sum(max(con_best-nonlcon_torlance,0),2);
    end
    if ~isempty(coneq_list)
        vio_best = vio_best+sum((abs(coneq_best)-nonlcon_torlance),2);
    end

    if INFORMATION_FLAG
        fprintf('%c  fval:    %f    violation:    %f    NFE:    %-3d\n',search_flag,fval_best,vio_best,NFE);
%         fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
%         if infor_search_flag == 0
%             fprintf('global x:          %s\n',num2str(x_global_infill));
%             fprintf('global value:      %f\n',fval_global_infill);
%             fprintf('global violation:  %s  %s\n',num2str(con_global_infill),num2str(coneq_global_infill));
%         else
%             fprintf('local  x:          %s\n',num2str(x_local_infill));
%             fprintf('local  value:      %f\n',fval_local_infill);
%             fprintf('local  violation:  %s  %s\n',num2str(con_local_infill),num2str(coneq_local_infill));
%         end
%         fprintf('\n');
    end
    
    result_x_best(iteration,:) = x_best;
    result_fval_best(iteration,:) = fval_best;
    iteration = iteration+1;
    
    % forced interrupt
    if iteration > iteration_max || NFE >= NFE_max
        done = 1;
    end
    
    % convergence judgment
    if CONVERGENCE_JUDGMENT_FLAG
        if (iteration > 2 && ...
                abs((fval_best-fval_best_old)/fval_best_old) < torlance)
            done = 1;
            if ~isempty(con_best)
                if sum(con_best > nonlcon_torlance)
                    done = 0;
                end
            end
            if ~isempty(coneq_best)
                if sum(abs(coneq_best) > nonlcon_torlance)
                    done = 0;
                end
            end
        end
    end

    fval_best_old = fval_best;
    search_flag = next_search_flag;
end

result_x_best = result_x_best(1:iteration-1,:);
result_fval_best = result_fval_best(1:iteration-1);

output.result_x_best = result_x_best;
output.result_fval_best = result_fval_best;


    function [fval,con,coneq] = modelFunction(x,object_function,nonlcon_function)
        % model function,concertrate fval,con,coneq into one function
        %
        if nargin < 3 || isempty(nonlcon_function)
            con = [];
            coneq = [];
        else
            [con,coneq] = nonlcon_function(x);
        end
        fval = object_function(x);
    end
    function [con,coneq] = totalconFunction...
            (x,nonlcon_function,cheapcon_function)
        con = [];
        coneq = [];
        if ~isempty(nonlcon_function)
            [expencon,expenconeq] = nonlcon_function(x);
            con = [con;expencon];
            coneq = [coneq;expenconeq];
        end
        if ~isempty(cheapcon_function)
            [expencon,expenconeq] = cheapcon_function(x);
            con = [con;expencon];
            coneq = [coneq;expenconeq];
        end
    end

    function X_new = differEvolutionRand(low_bou,up_bou,X,F,x_number,rand_number)
        if nargin < 4
            rand_number = 1;
            if nargin < 3
                x_number = 1;
                if nargin < 2
                    error('differEvolutionRand: lack scaling factor F');
                end
            end
        end
        [x_number__,variable_number__] = size(X);
        X_new = zeros(x_number,variable_number__);
        for x_index__ = 1:x_number
            index__ = randi(x_number__,2*rand_number+1,1);
            X_new(x_index__,:) = X(index__(1),:);
            for rand_index__ = 1:rand_number
                X_new(x_index__,:) = X_new(x_index__,:)+...
                    F*(X(index__(2*rand_index__),:)-X(index__(2*rand_index__+1),:));
                X_new(x_index__,:) = max(X_new(x_index__,:),low_bou);
                X_new(x_index__,:) = min(X_new(x_index__,:),up_bou);
            end
        end
    end
    function X_new = differEvolutionCurrentRand(low_bou,up_bou,X,F)
        [x_number__,variable_number__] = size(X);
        X_new = zeros(x_number__,variable_number__);
        for x_index__ = 1:x_number__
            index__ = randi(x_number__,3,1);
            X_new(x_index__,:) = X(x_index__,:)+...
                F*(X(index__(1),:)-X(x_index__,:)+...
                X(index__(2),:)-X(index__(3),:));
            X_new(x_index__,:) = max(X_new(x_index__,:),low_bou);
            X_new(x_index__,:) = min(X_new(x_index__,:),up_bou);
        end
    end
    function X_new = differEvolutionCurrentBest(low_bou,up_bou,X,F,x_best_index)
        [x_number__,variable_number__] = size(X);
        X_new = zeros(x_number__,variable_number__);
        for x_index__ = 1:x_number__
            index__ = randi(x_number__,2,1);
            X_new(x_index__,:) = X(x_index__,:)+...
                F*(X(x_best_index,:)-X(x_index__,:)+...
                X(index__(1),:)-X(index__(2),:));
            X_new(x_index__,:) = max(X_new(x_index__,:),low_bou);
            X_new(x_index__,:) = min(X_new(x_index__,:),up_bou);
        end
    end
    function X_new = differEvolutionCrossover(low_bou,up_bou,X,V,C_R)
        if size(X,1) ~= size(V,1)
            error('differEvolutionOffspring: size incorrect');
        end
        [x_number__,variable_number__] = size(X);
        X_new = X;
        rand_number = rand(x_number__,variable_number__);
        index__ = find(rand_number < C_R);
        X_new(index__) = V(index__);
        for x_index__ = 1:x_number__
            X_new(x_index__,:) = max(X_new(x_index__,:),low_bou);
            X_new(x_index__,:) = min(X_new(x_index__,:),up_bou);
        end
    end

    function [x_updata_list,fval_updata_list,con_updata_list,coneq_updata_list,NFE_updata,repeat_index] = dataLibraryUpdataProtect...
            (data_library_name,model_function,x_add_list,...
            x_list,low_bou,up_bou,protect_range)
        % function updata data with same_point_avoid protect
        % return fval
        % all list is x_number x variable_number matrix
        % notice if x_add is exist in library,point will be delete
        %
        variable_number__ = size(x_list,2);
        NFE_updata = 0;
        x_updata_list = [];fval_updata_list = [];con_updata_list = [];coneq_updata_list = [];repeat_index = [];
        for x_index__ = 1:size(x_add_list,1)
            x_updata__ = x_add_list(x_index__,:);
            
            % check x_potential if exist in data library
            % if not,updata data libraray
            distance__ = sum((abs(x_updata__-x_list)./(up_bou-low_bou)),2);
            [distance_min__,min_index__] = min(distance__);
            if distance_min__ < variable_number__*protect_range
                % distance to exist point of point to add is small than protect_range
                repeat_index = [repeat_index;min_index__];
            else
                [fval_updata__,con_updata__,coneq_updata__] = dataLibraryUpdata...
                    (data_library_name,model_function,x_updata__);NFE_updata = NFE_updata+1;
                x_updata_list = [x_updata_list;x_updata__];
                fval_updata_list = [fval_updata_list;fval_updata__];
                con_updata_list = [con_updata_list;con_updata__];
                coneq_updata_list = [coneq_updata_list;coneq_updata__];
            end
        end
    end

end

%% auxiliary function
function [x_best,fval_best,con_best,coneq_best] = findMinRaw...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance)
% find min fval in raw data
% x_list,rank is variable
% con_list,rank is con
% coneq_list,rank is coneq
% function will find min fval in con == 0
% if there was not feasible x,find min consum
%
con_best = [];
coneq_best = [];
max_nonlcon_list = zeros(size(x_list,1),1);
max_cheapcon_list = zeros(size(x_list,1),1);
% process expendsive con
if ~isempty(con_list)
    max_nonlcon_list = max(con_list,[],2);
end
if ~isempty(coneq_list)
    max_nonlcon_list = max(abs(coneq_list),[],2);
end

% add cheap con
for x_index = 1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq] = cheapcon_function(x_list(x_index,:));
        max_cheapcon_list(x_index) = max_cheapcon_list(x_index)+...
            sum(max(con,0))+sum(coneq.*coneq);
    end
end

con_judge_list = (max_nonlcon_list > nonlcon_torlance)+...
    (max_cheapcon_list > 0);
index = find(con_judge_list == 0);
if ~isempty(index)
    % feasible x
    x_list = x_list(index,:);
    fval_list = fval_list(index);
    if ~isempty(con_list)
        con_list = con_list(index,:);
    end
    if ~isempty(coneq_list)
        coneq_list = coneq_list(index,:);
    end

    % min fval
    [fval_best,index_best] = min(fval_list);
    x_best = x_list(index_best,:);
    if ~isempty(con_list)
        con_best = con_list(index_best,:);
    end
    if ~isempty(coneq_list)
        coneq_best = coneq_list(index_best,:);
    end
else
    % min consum
    [~,index_best] = min(max_nonlcon_list);
    fval_best = fval_list(index_best);
    x_best = x_list(index_best,:);
    if ~isempty(con_list)
        con_best = con_list(index_best,:);
    end
    if ~isempty(coneq_list)
        coneq_best = coneq_list(index_best,:);
    end
end
end

function [x_list,fval_list,con_list,coneq_list,vio_list,index_list] = rankData...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance)
% rank data base on feasibility rule
% infeasible is rank by sum of constraint
% torlance to cheapcon_function is 0
%
if nargin < 6 || isempty(nonlcon_torlance)
    nonlcon_torlance = 0;
end
if nargin < 5
    cheapcon_function = [];
end

[x_number,~] = size(x_list);
vio_list = zeros(x_number,1);
if ~isempty(con_list)
    vio_list = vio_list+sum(max(con_list-nonlcon_torlance,0),2);
end
if ~isempty(coneq_list)
    vio_list = vio_list+sum((abs(coneq_list)-nonlcon_torlance),2);
end

% add cheap con
for x_index = 1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq] = cheapcon_function(x_list(x_index,:));
        vio_list(x_index) = vio_list(x_index)+...
            sum(max(con,0))+sum(max(abs(coneq),0));
    end
end

% rank data
% infeasible data rank by violation, feasible data rank by fval
feasi_boolean_list = vio_list <= 0;
all = 1:x_number;
feasi_index_list = all(feasi_boolean_list);
infeasi_index_list = all(~feasi_boolean_list);
[~,index_list] = sort(fval_list(feasi_index_list));
feasi_index_list = feasi_index_list(index_list);
[~,index_list] = sort(vio_list(infeasi_index_list));
infeasi_index_list = infeasi_index_list(index_list);
index_list = [feasi_index_list,infeasi_index_list];

% rank by index_list
x_list = x_list(index_list,:);
fval_list = fval_list(index_list);
if ~isempty(con_list)
    con_list = con_list(index_list,:);
end
if ~isempty(coneq_list)
    coneq_list = coneq_list(index_list,:);
end
vio_list = vio_list(index_list);

end

%% surrogate model
function [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output] = getRadialBasisModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create radialbasis model and function
% if input model,function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con,coneq
% con is colume vector,coneq is colume vector
% var_function is same
%
basis_function = @(r) r.^3;
% basis_function = [];

[predict_function_fval,radialbasis_model_fval] = interpRadialBasisPreModel...
    (x_list,fval_list,basis_function);

if ~isempty(con_list)
    predict_function_con = cell(size(con_list,2),1);
    radialbasis_model_con = struct('X',[],'Y',[],...
        'beta',[],'radialbasis_matrix',[],'inv_radialbasis_matrix',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_con = repmat(radialbasis_model_con,[size(con_list,2),1]);
    for con_index = 1:size(con_list,2)
        [predict_function_con{con_index},radialbasis_model_con(con_index)] = interpRadialBasisPreModel...
            (x_list,con_list(:,con_index),basis_function);
    end
else
    predict_function_con = [];
    radialbasis_model_con = [];
end

if ~isempty(coneq_list)
    predict_function_coneq = cell(size(coneq_list,2),1);
    radialbasis_model_coneq = struct('X',[],'Y',[],...
        'beta',[],'radialbasis_matrix',[],'inv_radialbasis_matrix',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_coneq = repmat(radialbasis_model_coneq,[size(coneq_list,2),1]);
    for coneq_index = 1:size(coneq_list,2)
        [predict_function_coneq{coneq_index},radialbasis_model_con(con_index)] = interpRadialBasisPreModel...
            (x_list,coneq_list(:,coneq_index),basis_function);
    end
else
    predict_function_coneq = [];
    radialbasis_model_coneq = [];
end

object_function_surrogate = @(X_predict) objectFunctionSurrogate(X_predict,predict_function_fval);
if isempty(radialbasis_model_con) && isempty(radialbasis_model_coneq)
    nonlcon_function_surrogate = [];
else
    nonlcon_function_surrogate = @(X_predict) nonlconFunctionSurrogate(X_predict,predict_function_con,predict_function_coneq);
end

output.object_function_surrogate = object_function_surrogate;
output.nonlcon_function_surrogate = nonlcon_function_surrogate;

    function fval = objectFunctionSurrogate...
            (X_predict,predict_function_fval)
        % connect all predict favl
        %
        fval = predict_function_fval(X_predict);
    end
    function [con,coneq] = nonlconFunctionSurrogate...
            (X_predict,predict_function_con,predict_function_coneq)
        % connect all predict con and coneq
        %
        if isempty(predict_function_con)
            con = [];
        else
            con = zeros(size(X_predict,1),length(predict_function_con));
            for con_index__ = 1:length(predict_function_con)
                con(:,con_index__) = ....
                    predict_function_con{con_index__}(X_predict);
            end
        end
        if isempty(predict_function_coneq)
            coneq = [];
        else
            coneq = zeros(size(X_predict,1),length(predict_function_coneq));
            for coneq_index__ = 1:length(predict_function_coneq)
                coneq(:,coneq_index__) = ...
                    predict_function_coneq{coneq_index__}(X_predict);
            end
        end
    end
end

function [predict_function,radialbasis_model] = interpRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interp pre model function
% input initial data X,Y,which are real data
% X,Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model,include X,Y,base_function
% and predict_function
%
% Copyright 2023 Adel
%
if nargin < 3
    basis_function = [];
end

[x_number,variable_number] = size(X);

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__ = find(stdD_X == 0);
if ~isempty(index__),stdD_X(index__) = 1;end
index__ = find(stdD_Y == 0);
if ~isempty(index__),stdD_Y(index__) = 1;end
X_nomlz = (X-aver_X)./stdD_X;
Y_nomlz = (Y-aver_Y)./stdD_Y;

if isempty(basis_function)
    c = (prod(max(X_nomlz)-min(X_nomlz))/x_number)^(1/variable_number);
    basis_function = @(r) exp(-(r.^2)/c);
end

% initialization distance of all X
X_dis = zeros(x_number,x_number);
for variable_index = 1:variable_number
    X_dis = X_dis+(X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end
X_dis = sqrt(X_dis);

[beta,rdibas_matrix,inv_rdibas_matrix] = interpRadialBasis...
    (X_dis,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function = @(X_predict) interpRadialBasisPredictor...
    (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,beta,basis_function);

radialbasis_model.X = X;
radialbasis_model.Y = Y;
radialbasis_model.radialbasis_matrix = rdibas_matrix;
radialbasis_model.inv_radialbasis_matrix=inv_rdibas_matrix;
radialbasis_model.beta = beta;

radialbasis_model.aver_X = aver_X;
radialbasis_model.stdD_X = stdD_X;
radialbasis_model.aver_Y = aver_Y;
radialbasis_model.stdD_Y = stdD_Y;
radialbasis_model.basis_function = basis_function;

radialbasis_model.predict_function = predict_function;

% abbreviation:
% num: number,pred: predict,vari: variable
    function [beta,rdibas_matrix,inv_rdibas_matrix] = interpRadialBasis...
            (X_dis,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);
        
        % stabilize matrix
        rdibas_matrix = rdibas_matrix+eye(x_number)*1e-6;

        % get inverse matrix
        inv_rdibas_matrix = rdibas_matrix\eye(x_number);
        
        % solve beta
        beta = inv_rdibas_matrix*Y;
    end

    function [Y_pred] = interpRadialBasisPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,beta,basis_function)
        % radial basis function interpolation predict function
        %
        [x_pred_num,~] = size(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred-aver_X)./stdD_X;
        
        % calculate distance
        X_dis_pred = zeros(x_pred_num,x_num);
        for vari_index = 1:vari_num
            X_dis_pred = X_dis_pred+...
                (X_pred_nomlz(:,vari_index)-X_nomlz(:,vari_index)').^2;
        end
        X_dis_pred = sqrt(X_dis_pred);
        
        % predict variance
        Y_pred = basis_function(X_dis_pred)*beta;
        
        % normalize data
        Y_pred = Y_pred*stdD_Y+aver_Y;
    end

end

%% data library
function [fval_list,con_list,coneq_list] = dataLibraryUpdata...
    (data_library_name,model_function,x_list)
% updata data library
% updata format:
% variable_number,fval_number,con_number,coneq_number
% x,fval,con,coneq
%
fval_list = [];
con_list = [];
coneq_list = [];
[x_number,variable_number] = size(x_list);

if ~strcmp(data_library_name(end-3:end),'.txt')
    data_library_name = [data_library_name,'.txt'];
end

% store format
x_format_base = '%.8e ';
fval_format_base = '%.8e ';
x_format = repmat(x_format_base,1,variable_number);
file_optimalSurrogate_output = fopen(data_library_name,'a');
file_result = fopen('result_total.txt','a');

% updata format:
% variable_number,fval_number,con_number,coneq_number
% x,fval,con,coneq
for x_index = 1:x_number
    x = x_list(x_index,:);
    [fval,con,coneq] = model_function(x);
    fval = fval(:);
    con = con(:);
    coneq = coneq(:);
    fval_list = [fval_list;fval(:)'];
    con_list = [con_list;con(:)'];
    coneq_list = [coneq_list;coneq(:)'];
    
    % write data to txt_optimalSurrogateSADEKTS
    fprintf(file_optimalSurrogate_output,'%d ',variable_number);
    fprintf(file_optimalSurrogate_output,'%d ',length(fval));
    fprintf(file_optimalSurrogate_output,'%d ',length(con));
    fprintf(file_optimalSurrogate_output,'%d ',length(coneq));
    
    fprintf(file_optimalSurrogate_output,x_format,x);
    fval_format = repmat(fval_format_base,1,length(fval));
    fprintf(file_optimalSurrogate_output,fval_format,fval);
    fval_format = repmat(fval_format_base,1,length(con));
    fprintf(file_optimalSurrogate_output,fval_format,con);
    fval_format = repmat(fval_format_base,1,length(coneq));
    fprintf(file_optimalSurrogate_output,fval_format,coneq);
    fprintf(file_optimalSurrogate_output,'\n');
    
    % write data to txt_result
    fprintf(file_result,'%d ',variable_number);
    fprintf(file_result,'%d ',length(fval));
    fprintf(file_result,'%d ',length(con));
    fprintf(file_result,'%d ',length(coneq));
    
    fprintf(file_result,x_format,x);
    fval_format = repmat(fval_format_base,1,length(fval));
    fprintf(file_result,fval_format,fval);
    fval_format = repmat(fval_format_base,1,length(con));
    fprintf(file_result,fval_format,con);
    fval_format = repmat(fval_format_base,1,length(coneq));
    fprintf(file_result,fval_format,coneq);
    fprintf(file_result,'\n');
end

fclose(file_optimalSurrogate_output);
clear('file_optimalSurrogate_output');
fclose(file_result);
clear('file_result');
end

function [x_list,fval_list,con_list,coneq_list] = dataLibraryLoad...
    (data_library_name,low_bou,up_bou)
% load data from data library
% low_bou,up_bou is range of data
% updata format:
% variable_number,fval_number,con_number,coneq_number
% x,fval,con,coneq
%
if nargin < 3
    up_bou = inf;
    if nargin < 2
        low_bou = -inf;
        if nargin < 1
            error('dataLibraryLoad: lack data_library_name');
        end
    end
end

if ~strcmp(data_library_name(end-3:end),'.txt')
    data_library_name = [data_library_name,'.txt'];
end

% updata format:
% variable_number,fval_number,con_number,coneq_number
% x,fval,con,coneq
if exist(data_library_name,'file') == 2
    data_list = importdata(data_library_name);
    if ~isempty(data_list)
        % search whether exist point
        x_list = [];
        fval_list = [];
        con_list = [];
        coneq_list = [];
        
        for data_index = 1:size(data_list,1)
            data = data_list(data_index,:);
            
            variable_number = data(1);
            fval_number = data(2);
            con_number = data(3);
            coneq_number = data(4);
            
            base = 5;
            x = data(base:base+variable_number-1);
            judge = sum(x < low_bou)+sum(x > up_bou);
            if ~judge
                x_list = [x_list;x];
                base = base+variable_number;
                fval_list = [fval_list;data(base:base+fval_number-1)];
                base = base+fval_number;
                con = data(base:base+con_number-1);
                if ~isempty(con)
                    con_list = [con_list;con];
                end
                base = base+con_number;
                coneq = data(base:base+coneq_number-1);
                if ~isempty(coneq)
                    coneq_list = [coneq_list;coneq];
                end
            end
        end
    else
        x_list = [];
        fval_list = [];
        con_list = [];
        coneq_list = [];
    end
else
    x_list = [];
    fval_list = [];
    con_list = [];
    coneq_list = [];
end
end

%% LHD
function [X,X_new,distance_min_nomlz] = getLatinHypercube...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou,cheapcon_function)
% generate sample sequence latin hypercube
% election sequential method is used(sample and iteration)
%
% sample number is total point in area
% default low_bou is 0,up_bou is 1,cheapcon_function is []
% low_bou and up_bou is colume vector
% X_exist,X,supply_X is x_number x variable_number matrix
% X_exist should meet bou
%
% reference:
% [1]LONG T,LI X,SHI R,et al.,Gradient-Free Trust-Region-Based
% Adaptive Response Surface Method for Expensive Aircraft Optimization[J].
% AIAA Journal,2018,56(2): 862-73.
%
% Copyright 2022 Adel
%
if nargin < 6
    cheapcon_function = [];
    if nargin < 5
        if nargin < 3
            if nargin < 2
                error('getLatinHypercube: lack variable_number');
            end
        end
        low_bou = zeros(variable_number,1);
        up_bou = ones(variable_number,1);
    end
end

% check x_exist_list if meet boundary
if ~isempty(X_exist)
    if size(X_exist,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
    index = find(X_exist < low_bou);
    index = [index,find(X_exist > up_bou)];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list range error');
    end
    X_exist_nomlz = (X_exist-low_bou)./(up_bou-low_bou);
else
    X_exist_nomlz = [];
end

if sample_number <= 0
    X = [];
    X_new = [];
    distance_min_nomlz = [];
    return;
end

iteration_max = 1000*variable_number;
x_new_number = sample_number-size(X_exist,1);
if x_new_number <= 0
    X = X_exist;
    X_new = [];
    distance_min_nomlz = getMinDistance(X_exist_nomlz);
    return;
end

% get quasi-feasible point
x_initial_number = 100*x_new_number;
if ~isempty(cheapcon_function)
    X_supply_quasi_nomlz = [];
    % check if have enough X_supply_nomlz
    while size(X_supply_quasi_nomlz,1) < 100*x_new_number
        X_supply_initial_nomlz = rand(x_initial_number,variable_number);
        x_index = 1;
        while x_index <= size(X_supply_initial_nomlz,1)
            x_supply = X_supply_initial_nomlz(x_index,:).*(up_bou-low_bou)+low_bou;
            if cheapcon_function(x_supply) > 0
                X_supply_initial_nomlz(x_index,:) = [];
            else
                x_index = x_index+1;
            end
        end
        X_supply_quasi_nomlz = [X_supply_quasi_nomlz;X_supply_initial_nomlz];
    end
else
    X_supply_quasi_nomlz = rand(x_initial_number,variable_number);
end

% iterate and get final x_supply_list
iteration = 0;
x_supply_quasi_number = size(X_supply_quasi_nomlz,1);
distance_min_nomlz = 0;
X_new_nomlz = [];
while iteration <= iteration_max
    % random select x_new_number X to X_trial_nomlz
    x_select_index = randperm(x_supply_quasi_number,x_new_number);
    
    % get distance min itertion X_
    distance_min_iteration = getMinDistanceIter...
        (X_supply_quasi_nomlz(x_select_index,:),X_exist_nomlz);
    
    % if distance_min_iteration is large than last time
    if distance_min_iteration > distance_min_nomlz
        distance_min_nomlz = distance_min_iteration;
        X_new_nomlz = X_supply_quasi_nomlz(x_select_index,:);
    end
    
    iteration = iteration+1;
end
X_new = X_new_nomlz.*(up_bou-low_bou)+low_bou;
X = [X_new;X_exist];
distance_min_nomlz = getMinDistance([X_new_nomlz;X_exist_nomlz]);

    function distance_min__ = getMinDistance(x_list__)
        % get distance min from x_list
        %
        % sort x_supply_list_initial to decrese distance calculate times
        x_list__ = sortrows(x_list__,1);
        [sample_number__,variable_number__] = size(x_list__);
        distance_min__ = variable_number__;
        for x_index__ = 1:sample_number__
            x_curr__ = x_list__(x_index__,:);
            x_next_index__ = x_index__ + 1;
            % only search in min_distance(x_list had been sort)
            search_range__ = variable_number__;
            while x_next_index__ <= sample_number__ &&...
                    (x_list__(x_next_index__,1)-x_list__(x_index__,1))^2 ...
                    < search_range__
                x_next__ = x_list__(x_next_index__,:);
                distance_temp__ = sum((x_next__-x_curr__).^2);
                if distance_temp__ < distance_min__
                    distance_min__ = distance_temp__;
                end
                if distance_temp__ < search_range__
                    search_range__ = distance_temp__;
                end
                x_next_index__ = x_next_index__+1;
            end
        end
        distance_min__ = sqrt(distance_min__);
    end
    function distance_min__ = getMinDistanceIter...
            (x_list__,x_exist_list__)
        % get distance min from x_list
        %
        % sort x_supply_list_initial to decrese distance calculate times
        x_list__ = sortrows(x_list__,1);
        [sample_number__,variable_number__] = size(x_list__);
        distance_min__ = variable_number__;
        for x_index__ = 1:sample_number__
            x_curr__ = x_list__(x_index__,:);
            x_next_index__ = x_index__ + 1;
            % only search in min_distance(x_list had been sort)
            search_range__ = variable_number__;
            while x_next_index__ <= sample_number__ &&...
                    (x_list__(x_next_index__,1)-x_list__(x_index__,1))^2 ...
                    < search_range__
                x_next__ = x_list__(x_next_index__,:);
                distance_temp__ = sum((x_next__-x_curr__).^2);
                if distance_temp__ < distance_min__
                    distance_min__ = distance_temp__;
                end
                if distance_temp__ < search_range__
                    search_range__ = distance_temp__;
                end
                x_next_index__ = x_next_index__+1;
            end
            for x_exist_index = 1:size(x_exist_list__,1)
                x_next__ = x_exist_list__(x_exist_index,:);
                distance_temp__ = sum((x_next__-x_curr__).^2);
                if distance_temp__ < distance_min__
                    distance_min__ = distance_temp__;
                end
            end
        end
        distance_min__ = sqrt(distance_min__);
    end
end
