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

% variable_number = 20;
% object_function = @(x) benchmark.singleEP20Object(x);
% object_function_LF = @(x) benchmark.singleEP20ObjectLow(x);
% A = [];
% B = [];
% Aeq = [];
% Beq = [];
% low_bou = ones(1,variable_number)*-30;
% up_bou = ones(1,variable_number)*30;
% nonlcon_function = [];
% nonlcon_function_LF = [];
% cheapcon_function = [];

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

variable_number = 2;
object_function = @(x) benchmark.singleG06Object(x);
object_function_LF = @(x) benchmark.singleG06ObjectLow(x);
A = [];
B = [];
Aeq = [];
Beq = [];
low_bou = [13,0];
up_bou = [100,100];
nonlcon_function = @(x) benchmark.singleG06Nonlcon(x);
nonlcon_function_LF = @(x) benchmark.singleG06NonlconLow(x);
cheapcon_function = [];
model_function = @(x) modelFunction(x,@(x) benchmark.singleG01Object(x),@(x) benchmark.singleG06Nonlcon(x));
model_function_LF = @(x) modelFunction(x,@(x) benchmark.singleG01Object(x),@(x) benchmark.singleG06Nonlcon(x));

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
% [x_best,fval_best,~,output] = fmincon(object_function,x_initial,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',10000,'Display','iter-detailed'))

%% single run

delete([data_library_name,'.txt']);
delete('result_total.txt');

LF_HF_ratio = 3;

[x_best,fval_best,NFE,output] = optimalSurrogateMFRBFMFGPC...
    (model_function,variable_number,low_bou,up_bou,...
    cheapcon_function,model_function_LF,LF_HF_ratio,30,500)
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
function [x_best,fval_best,NFE,output] = optimalSurrogateMFRBFMFGPC...
    (model_function,model_function_LF,LF_HF_ratio,variable_number,low_bou,up_bou,...
    cheapcon_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance,x_initial_list)
% MFRBF-MFGPC optimization algorithm
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
    cheapcon_function = [];
end

DRAW_FIGURE_FLAG = 0; % whether draw data
INFORMATION_FLAG = 1; % whether print data
CONVERGENCE_JUDGMENT_FLAG = 1; % whether judgment convergence

% hyper parameter
sample_number_initial = 2+variable_number;
sample_number_initial_LF = LF_HF_ratio*sample_number_initial;

trial_number = min(100*variable_number,100);
coord_select_prob_initial = min(20/variable_number,1);
sigma_coord_initial = 0.1*(up_bou-low_bou);
sigma_coord_min = 0.2*1/64*(up_bou-low_bou);
sigma_coord_max = 2*(up_bou-low_bou);
tau_success = 3;
tau_fail = max(variable_number,5);

% max fval when normalize fval,con,coneq
nomlz_fval = 10;

protect_range = 1e-8;

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
%     if expensive_nonlcon_flag
%         NFE_max = 50*variable_number;
%     else
%         NFE_max = 20*variable_number;
%     end
    NFE_max = 10+10*variable_number;
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
C_repeat = 0;
next_search_flag = 'G'; % 'G' is global search,'l' is local search
hyp.mean = 0;
hyp.cov = [0,0];
while ~done
    search_flag = next_search_flag;
    % nomalization all data by max fval
    fval_max = max(abs(fval_list),[],1);
    fval_nomlz_list = fval_list/fval_max*nomlz_fval;
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
    [RBF_model_fval,RBF_model_con,RBF_model_coneq,output_RBF] = getMFRBFModel...
        (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list);
    object_function_surrogate = output_RBF.object_function_surrogate;
    nonlcon_function_surrogate = output_RBF.nonlcon_function_surrogate;

    %     R2 = 1-sum((RBF_model_fval.beta./diag(RBF_model_fval.inv_radialbasis_matrix)).^2)/sum((RBF_model_fval.Y-mean(RBF_model_fval.Y)).^2);
    %     fprintf('global R2:  %f\n',R2);

    improve_flag = 0;
%     search_flag = 'G';
    if search_flag == 'G'
        % global search

        % calculate gradient
        [~,gradient] = differ(object_function_surrogate,x_best,[],1e-3);
        gradient = abs(gradient);
        gradient_max = max(gradient); gradient_min = min(gradient);
        gradient = (gradient-gradient_min)./(gradient_max-gradient_min);

        % step 5
        % perturbed to get trial point
        coord_select_prob = coord_select_prob_initial*(1-log(NFE-sample_number_initial+1)/log(NFE_max-sample_number_initial));

        % select coordinates to perturb
        coordinate_select_list = [];
        for variable_index = 1:variable_number
            if rand() < coord_select_prob*gradient(variable_index)
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
            %             % using KS equation to assemble constraints
            %             KS_nomlz_list = max(max(con_list,[],2),max(coneq_list,[],2));
            %
            %             [RBF_model_KS,~,~,output_RBF] = getRadialBasisModel...
            %                 (x_list,KS_nomlz_list,[],[]);
            %             KS_function_surrogate = output_RBF.object_function_surrogate;
            %
            %             [KS_pred_trial_list] = KS_function_surrogate(x_trial_list);
            %             KS_var_trial_list = fval_var_trial_list;

            %             vio_trial_list = zeros(trial_number,1);
            %             if ~isempty(con_nomlz_list)
            %                 vio_trial_list = vio_trial_list+sum(max(con_pred_trial_list-nonlcon_torlance,0),2);
            %             end
            %             if ~isempty(coneq_nomlz_list)
            %                 vio_trial_list = vio_trial_list+sum((abs(con_pred_trial_list)-nonlcon_torlance),2);
            %             end

            %             feasi_boolean_trial_list = vio_trial_list <= nonlcon_torlance;
            %
            %             % if have feasiable_index_list,only use feasiable to choose
            %             if all(~feasi_boolean_trial_list)
            %                 % base on constaints improve select global infill
            %                 % lack process of equal constraints
            %                 con_nomlz_base = max(min(con_nomlz_list,[],1),0);
            %                 con_impove_probability_list = sum(...
            %                     normcdf((con_nomlz_base-con_pred_trial_list)./sqrt(con_var_trial_list)),2);
            %                 [~,con_best_index] = max(con_impove_probability_list);
            %                 con_best_index = con_best_index(1);
            %                 x_global_infill = x_trial_list(con_best_index,:);
            %             else
            %                 % base on fitness DE point to select global infill
            %                 if expensive_nonlcon_flag
            %                     x_trial_list = x_trial_list(feasi_boolean_trial_list,:);
            %                     fval_pred_trial_list = fval_pred_trial_list(feasi_boolean_trial_list);
            %                     fval_var_trial_list = fval_var_trial_list(feasi_boolean_trial_list);
            %                 end
            %             end

            % base on filter to decide which x should be choose
            pareto_index_list = getParetoFront([fval_list,vio_list]);
            class_list = -1*ones(size(x_list,1),1);
            class_list(pareto_index_list) = 1;

            [predict_function,CGP_model] = classifyGaussProcess(x_list,class_list,hyp);
            [~,pro_pred] = predict_function(x_trial_list);

            fval_pred_trial_min = min(fval_pred_trial_list,[],1);
            fval_pred_trial_max = max(fval_pred_trial_list,[],1);
            fval_var_trial_min = min(fval_var_trial_list,[],1);
            fval_var_trial_max = max(fval_var_trial_list,[],1);

            fitness_list = -(fval_pred_trial_list-fval_pred_trial_min)/(fval_pred_trial_max-fval_pred_trial_min)-...
                (fval_var_trial_max-fval_var_trial_list)/(fval_var_trial_max-fval_var_trial_min);
            fitness_list = fitness_list.*pro_pred;

            [~,fitness_best_index] = min(fitness_list);
            fitness_best_index = fitness_best_index(1);
            x_global_infill = x_trial_list(fitness_best_index,:);

        else
            fval_pred_trial_min = min(fval_pred_trial_list,[],1);
            fval_pred_trial_max = max(fval_pred_trial_list,[],1);
            fval_var_trial_min = min(fval_var_trial_list,[],1);
            fval_var_trial_max = max(fval_var_trial_list,[],1);

            fitness_list = -(fval_pred_trial_list-fval_pred_trial_min)/(fval_pred_trial_max-fval_pred_trial_min)-...
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
        vio_global_infill = 0;
        if ~isempty(con_list)
            con_list = [con_list;con_global_infill];
            vio_global_infill = vio_global_infill+sum(max(con_global_infill-nonlcon_torlance,0),2);
        end
        if ~isempty(coneq_list)
            coneq_list = [coneq_list;coneq_global_infill];
            vio_global_infill = vio_global_infill+sum((abs(coneq_global_infill)-nonlcon_torlance),2);
        end
        vio_list = [vio_list;vio_global_infill];

        % whether impove pupolation,if imporve,continue global
        % notice last one is x_local fval,con and coneq
        next_search_flag = 'l';
        if isempty(x_global_infill)
            x_global_infill = x_list(repeat_index,:);
            fval_global_infill = fval_list(repeat_index,:);
            if ~isempty(con_list)
                con_global_infill = con_list(repeat_index,:);
            end
            if ~isempty(coneq_list)
                coneq_global_infill = coneq_list(repeat_index,:);
            end
            vio_global_infill = vio_list(repeat_index,:);
        else
            % judge if improve solve
            if all(~feasi_boolean_list)
                min_vio = min(vio_list(1:end-1));

                % improve,continue global search
                if vio_global_infill < min_vio
                    next_search_flag = 'G';
                    improve_flag = 1;
                end
            else
                if expensive_nonlcon_flag
                    min_fval = min(fval_list([feasi_boolean_list;false(1)]),[],1);
                    if vio_global_infill <= 0 % if exist constraint,infill point should feasible
                        vio_judge = 1;
                    else
                        vio_judge = 0;
                    end
                else
                    min_fval = min(fval_list(1:end-1));
                    vio_judge = 1;
                end

                % imporve,continue global search
                if fval_global_infill < min_fval && vio_judge
                    next_search_flag = 'G';
                    improve_flag = 1;
                end
            end

            if DRAW_FIGURE_FLAG && variable_number < 3
                interpVisualize(RBF_model_fval,low_bou,up_bou);
                line(x_global_infill(1),x_global_infill(2),fval_global_infill./fval_max*nomlz_fval,'Marker','o','color','r');
                hold on;
                scatter3(x_trial_list(:,1),x_trial_list(:,2),fitness_list*10,'black');
                hold off;
            end

        end

    else
        % local search
        [~,~,~,~,~,index_list] = rankData...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            cheapcon_function,nonlcon_torlance);

        % step 8
        % rand select initial local point from x_list

        x_index = index_list(randi(min(size(x_list,1),2*variable_number)));
%         x_index = index_list(1);
        x_initial = x_list(x_index,:);

        %         % select nearest point to construct RBF
        %         distance = sum(((x_initial-x_list)./(up_bou-low_bou)).^2,2);
        %         [~,index_list] = sort(distance);
        %
        %         x_number = size(x_list,1);
        %         R2_list = zeros(x_number,1);
        %         for x_index = (2*variable_number):x_number
        %             RBF_number = x_index;
        %
        %             index_local_list = index_list(1:RBF_number);
        %             X = x_list(index_local_list,:);
        %             Y = fval_nomlz_list(index_local_list,:);
        %
        %             basis_function = @(r) r.^3;
        %             [predict_function,radialbasis_model] = interpRadialBasisPreModel(X,Y,basis_function);
        %             error_list = radialbasis_model.beta./diag(radialbasis_model.inv_radialbasis_matrix);
        %
        %             R2_list(x_index) = 1-sum(error_list.^2)/sum((Y-mean(Y)).^2);
        %
        %             if R2_list(x_index) > 0.95
        %                 break;
        %             end
        %         end
        %         fprintf('local number: %d \n',x_index);
        %
        %         x_local_list = x_list(index_local_list,:);
        %         fval_local_nomlz_list = fval_nomlz_list(index_local_list,:);
        %         if ~isempty(con_list)
        %             con_local_nomlz_list = con_nomlz_list(index_local_list,:);
        %         else
        %             con_local_nomlz_list = [];
        %         end
        %         if ~isempty(coneq_list)
        %             coneq_local_nomlz_list = coneq_nomlz_list(index_local_list,:);
        %         else
        %             coneq_local_nomlz_list = [];
        %         end
        %
        %         [local_model_fval,local_model_con,local_model_coneq,output_local] = getRadialBasisModel...
        %             (x_local_list,fval_local_nomlz_list,con_local_nomlz_list,coneq_local_nomlz_list);
        %         object_function_surrogate = output_local.object_function_surrogate;
        %         nonlcon_function_surrogate = output_local.nonlcon_function_surrogate;

%         [RBF_model_fval,RBF_model_con,RBF_model_coneq,output_RBF] = getEnsemleRadialBasisModel...
%             (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list);
%         object_function_surrogate = output_RBF.object_function_surrogate;
%         nonlcon_function_surrogate = output_RBF.nonlcon_function_surrogate;

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

        x_list = [x_list;x_local_infill];
        fval_list = [fval_list;fval_local_infill];
        vio_local_infill = 0;
        if ~isempty(con_list)
            con_list = [con_list;con_local_infill];
            vio_local_infill = vio_local_infill+sum(max(con_local_infill-nonlcon_torlance,0),2);
        end
        if ~isempty(coneq_list)
            coneq_list = [coneq_list;coneq_local_infill];
            vio_local_infill = vio_local_infill+sum((abs(coneq_local_infill)-nonlcon_torlance),2);
        end
        vio_list = [vio_list;vio_local_infill];

        % whether impove pupolation,if imporve,continue local
        % notice last one is x_local fval,con and coneq
        next_search_flag = 'G';
        if isempty(x_local_infill)
            x_local_infill = x_list(repeat_index,:);
            fval_local_infill = fval_list(repeat_index,:);
            if ~isempty(con_list)
                con_local_infill = con_list(repeat_index,:);
            end
            if ~isempty(coneq_list)
                coneq_local_infill = coneq_list(repeat_index,:);
            end
            vio_local_infill = vio_list(repeat_index,:);
        else
            if all(~feasi_boolean_list)
                min_vio = min(vio_list(1:end-1));

                % improve,continue local search
                if vio_local_infill < min_vio
                    next_search_flag = 'l';
                    improve_flag = 1;
                end
            else
                if expensive_nonlcon_flag
                    min_fval = min(fval_list([feasi_boolean_list;false(1)]),[],1);
                    if vio_local_infill <= 0 % if exist constraint,infill point should feasible
                        vio_judge = 1;
                    else
                        vio_judge = 0;
                    end
                else
                    min_fval = min(fval_list(1:end-1));
                    vio_judge = 1;
                end

                % imporve,continue local search
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

    end

    if search_flag == 'G'
        fval_updata = fval_global_infill;
        con_updata = con_global_infill;
        coneq_updata = coneq_global_infill;
        vio_updata = vio_global_infill;
    else
        fval_updata = fval_local_infill;
        con_updata = con_local_infill;
        coneq_updata = coneq_local_infill;
        vio_updata = vio_local_infill;
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
%         if search_flag == 'G'
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
                abs((fval_updata-fval_updata_old)/fval_updata_old) < torlance)
            if vio_updata == 0
                C_repeat = C_repeat+1;
                if C_repeat == 2
                    %                     done = 1;

%                     x_updata_list = lhsdesign(sample_number_initial,variable_number).*(up_bou-low_bou)+low_bou;
%                     [fval_updata_list,con_updata_list,coneq_updata_list] = dataLibraryUpdata...
%                         (data_library_name,model_function,x_updata_list);NFE = NFE+size(x_updata_list,1);
                end
            end
        else
            C_repeat = 0;
        end
    end

    fval_best_old = fval_best;

    fval_updata_old = fval_updata;
    con_updata_old = con_updata;
    coneq_updata_old = coneq_updata;
    vio_updata_old = vio_updata;
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
% infeasible data rank by violation,feasible data rank by fval
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

function pareto_index_list = getParetoFront(data_list)
% distinguish pareto front of data list
% data_list is x_number x data_number matrix
% notice if all data of x1 is less than x2,x1 domain x2
%
x_number = size(data_list,1);
pareto_index_list = []; % sort all index of filter point list

% select no domain filter
for x_index = 1:x_number
    data = data_list(x_index,:);
    pareto_index = 1;
    add_filter_flag = 1;
    while pareto_index <= length(pareto_index_list)
        % compare x with exit pareto front point
        x_pareto_index = pareto_index_list(pareto_index,:);

        % contain constraint of x_filter
        data_pareto = data_list(x_pareto_index,:);

        % compare x with x_pareto
        judge = data >= data_pareto;
        if ~sum(~judge)
            add_filter_flag = 0;
            break;
        end

        % if better or equal than exit pareto point,reject pareto point
        judge = data <= data_pareto;
        if ~sum(~judge)
            pareto_index_list(pareto_index) = [];
            pareto_index = pareto_index-1;
        end

        pareto_index = pareto_index+1;
    end

    % add into pareto list if possible
    if add_filter_flag
        pareto_index_list = [pareto_index_list;x_index];
    end
end
end

%% GPC
function [predict_function,CGPMF_model] = classifyGaussProcessMultiFidelity...
    (XHF,ClassHF,XLF,ClassLF,hyp)
% generate multi fidelity gaussian process classifier model
% version 6,this version is assembly of gpml-3.6 EP method
% X,XHF,XLF is x_number x variable_number matirx
% Class,ClassH,ClassLF is x_number x 1 matrix
% low_bou,up_bou is 1 x variable_number matrix
%
% input:
% XHF,ClassHF,XLF,ClassLF,hyp(mean,cov(lenD,etaD,lenL,etaL,rho))
%
% abbreviation:
% pred: predicted,nomlz: normalization,num: number
% var: variance
%
% reference: [1]FSCAB C,PP D,EK E,et al. Multi-fidelity classification
% using Gaussian processes: Accelerating the prediction of large-scale
% computational models [J]. Computer Methods in Applied Mechanics and
% Engineering,357(C): 112602-.
%
X = [XHF;XLF];
Class = [ClassHF;ClassLF];
[x_number,variable_number] = size(X);
x_HF_number = size(XHF,1);
x_LF_number = size(XLF,1);
if nargin < 5
    hyp.mean = 0;
    hyp.cov = zeros(1,variable_number*2+3);
end

% normalization data
aver_X = mean(X);
stdD_X = std(X);
index__ = find(stdD_X  ==  0);
if  ~isempty(index__),stdD_X(index__) = 1; end
X_nomlz = (X-aver_X)./stdD_X;

object_function = @(x) objectNLLMFGPC(x,{@infEP},{@meanConst},{@calCovMF},{@likErf},X_nomlz,Class);
hyp_x = [hyp.mean,hyp.cov];

% [fval,gradient] = object_function(hyp_x)
% [fval_differ,gradient_differ] = differ(object_function,hyp_x)

low_bou_hyp = -4*ones(1,2*variable_number+4);
up_bou_hyp = 4*ones(1,2*variable_number+4);
hyp_x = fmincon(object_function,hyp_x,[],[],[],[],low_bou_hyp,up_bou_hyp,[],...
    optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',true,...
    'MaxFunctionEvaluations',20,'OptimalityTolerance',1e-6));

hyp.mean = hyp_x(1);
hyp.cov = hyp_x(2:end);
hyp.lik = [];
post = infEP(hyp,{@meanConst},{@calCovMF},{@likErf},X_nomlz,Class);
predict_function = @(x_pred) classifyGaussPredictor...
    (x_pred,hyp,{@meanConst},{@calCovMF},{@likErf},post,X_nomlz,aver_X,stdD_X);

% output model
X = {XHF,XLF};
Class = {ClassHF,ClassLF};
CGPMF_model.X = X;
CGPMF_model.Class = Class;
CGPMF_model.aver_X = aver_X;
CGPMF_model.stdD_X = stdD_X;
CGPMF_model.predict_function = predict_function;
CGPMF_model.hyp = hyp;
CGPMF_model.post = post;

    function [fval,gradient] = objectNLLMFGPC(x,inf,mean,cov,lik,X,Y)
        hyp_iter.mean = x(1);
        hyp_iter.cov = x(2:end);
        hyp_iter.lik = [];

        if nargout < 2
            [~,nlZ] = feval(inf{:},hyp_iter,mean,cov,lik,X,Y);
            fval = nlZ;
        elseif nargout < 3
            [~,nlZ,dnlZ] = feval(inf{:},hyp_iter,mean,cov,lik,X,Y);
            fval = nlZ;
            gradient = [dnlZ.mean,dnlZ.cov];
        end
    end

    function [class,possibility,miu_pre,var_pre] = classifyGaussPredictor...
            (X_pred,hyp,mean,cov,lik,post,X,aver_X,stdD_X)
        % predict function
        %
        X_pred_nomlz = (X_pred-aver_X)./stdD_X;
        pred_num = size(X_pred_nomlz,1);
        ys = ones(pred_num,1);

        alpha = post.alpha; L = post.L; sW = post.sW;
        %verify whether L contains valid Cholesky decomposition or something different
        Lchol = isnumeric(L) && all(all(tril(L,-1) == 0)&diag(L)'>0&isreal(diag(L))');
        ns = size(X_pred_nomlz,1);                                       % number of data points
        nperbatch = 1000;                       % number of data points per mini batch
        nact = 0;                       % number of already processed test data points
        ymu = zeros(ns,1); ys2 = ymu; miu_pre = ymu; var_pre = ymu; possibility = ymu;   % allocate mem
        while nact<ns               % process minibatches of test cases to save memory
            id = (nact+1):min(nact+nperbatch,ns);               % data points to process
            kss = feval(cov{:},hyp.cov,X_pred_nomlz(id,:),'diag');              % self-variance
            Ks = feval(cov{:},hyp.cov,X,X_pred_nomlz(id,:));        % avoid computation
            ms = feval(mean{:},hyp.mean,X_pred_nomlz(id,:));
            N = size(alpha,2);  % number of alphas (usually 1; more in case of sampling)
            Fmu = repmat(ms,1,N) + Ks'*full(alpha);        % conditional mean fs|f
            miu_pre(id) = sum(Fmu,2)/N;                                   % predictive means
            if Lchol    % L contains chol decomp  = > use Cholesky parameters (alpha,sW,L)
                V  = L'\(repmat(sW,1,length(id)).*Ks);
                var_pre(id) = kss - sum(V.*V,1)';                       % predictive variances
            else                % L is not triangular  = > use alternative parametrisation
                if isnumeric(L),LKs = L*Ks; else LKs = L(Ks); end    % matrix or callback
                var_pre(id) = kss + sum(Ks.*LKs,1)';                    % predictive variances
            end
            var_pre(id) = max(var_pre(id),0);   % remove numerical noise i.e. negative variances
            Fs2 = repmat(var_pre(id),1,N);     % we have multiple values in case of sampling
            if nargin<9
                [Lp,Ymu,Ys2] = feval(lik{:},hyp.lik,[],Fmu(:),Fs2(:));
            else
                Ys = repmat(ys(id),1,N);
                [Lp,Ymu,Ys2] = feval(lik{:},hyp.lik,Ys(:),Fmu(:),Fs2(:));
            end
            possibility(id)  = sum(reshape(Lp,[],N),2)/N;    % log probability; sample averaging
            ymu(id) = sum(reshape(Ymu,[],N),2)/N;          % predictive mean ys|y and ..
            ys2(id) = sum(reshape(Ys2,[],N),2)/N;                          % .. variance
            nact = id(end);          % set counter to index of last processed data point
        end

        possibility = exp(possibility);
        class = ones(pred_num,1);
        index_list = possibility < 0.5;
        class(index_list) = -1;
    end

    function [K,dK_dcov] = calCovMF(cov,X,Z)
        % obtain covariance of x
        % cov: lenD,etaD,lenL,etaL,rho
        % len equal to 1/len_origin.^2
        %
        % % k = eta*exp(-sum(x_dis*theta)/vari_num);
        %
        [x_num,vari_num] = size(X);

        lenD = exp(cov(1:vari_num));
        etaD = exp(cov(vari_num+1));
        lenL = exp(cov(vari_num+1+(1:vari_num)));
        etaL = exp(cov(2*(vari_num+1)));
        rho = exp(cov(end));

        if nargin > 2 && nargout < 2 && ~isempty(Z)
            if strcmp(Z,'diag')
                K = rho*rho*etaL+etaD;
                return
            end
        end

        % predict
        if nargin > 2 && nargout < 2 && ~isempty(Z)
            [z_num,vari_num] = size(Z);
            % initializate square of X inner distance sq/ vari_num
            sq_dis_v = zeros(x_num,z_num,vari_num);
            for len_index = 1:vari_num
                sq_dis_v(:,:,len_index) = (X(:,len_index)-Z(:,len_index)').^2/vari_num;
            end

            % exp of x__x with D
            exp_disD = zeros(x_HF_number,z_num);
            for len_index = 1:vari_num
                exp_disD = exp_disD+...
                    sq_dis_v(1:x_HF_number,:,len_index)*lenD(len_index);
            end
            exp_disD = exp(-exp_disD);

            % exp of x__x with L
            exp_disL = zeros(x_num,z_num);
            for len_index = 1:vari_num
                exp_disL = exp_disL+...
                    sq_dis_v(1:x_num,:,len_index)*lenL(len_index);
            end
            exp_disL = exp(-exp_disL);

            K = exp_disL;
            K(1:x_HF_number,:) = rho*rho*etaL*K(1:x_HF_number,:)+etaD*exp_disD;
            K(x_HF_number+1:end,:) = rho*etaL*K(x_HF_number+1:end,:);
        else
            % initializate square of X inner distance sq/ vari_num
            sq_dis_v = zeros(x_num,x_num,vari_num);
            for len_index = 1:vari_num
                sq_dis_v(:,:,len_index) = (X(:,len_index)-X(:,len_index)').^2/vari_num;
            end

            % exp of x__x with H
            exp_disD = zeros(x_num);
            for len_index = 1:vari_num
                exp_disD(1:x_HF_number,1:x_HF_number) = exp_disD(1:x_HF_number,1:x_HF_number)+...
                    sq_dis_v(1:x_HF_number,1:x_HF_number,len_index)*lenD(len_index);
            end
            exp_disD(1:x_HF_number,1:x_HF_number) = exp(-exp_disD(1:x_HF_number,1:x_HF_number));
            KD = etaD*exp_disD;

            % exp of x__x with L
            exp_disL = zeros(x_num);
            for len_index = 1:vari_num
                exp_disL = exp_disL+...
                    sq_dis_v(1:end,1:end,len_index)*lenL(len_index);
            end
            exp_disL = exp(-exp_disL);
            eta_exp_disL = etaL*exp_disL;

            % times rho: HH to rho2,HL to rho,LL to 1
            KL = eta_exp_disL;
            KL(1:x_HF_number,1:x_HF_number) = ...
                (rho*rho)*eta_exp_disL(1:x_HF_number,1:x_HF_number);
            KL(1:x_HF_number,(x_HF_number+1):end) = ...
                rho*eta_exp_disL(1:x_HF_number,(x_HF_number+1):end);
            KL((x_HF_number+1):end,1:x_HF_number) = ...
                KL(1:x_HF_number,(x_HF_number+1):end)';

            K = KL+KD;

            if nargout >= 2
                dK_dcov = cell(1,2*vari_num+3);

                % len D
                for len_index = 1:vari_num
                    dK_dlenD = zeros(x_num);
                    dK_dlenD(1:x_HF_number,1:x_HF_number) = -KD(1:x_HF_number,1:x_HF_number).*...
                        sq_dis_v(1:x_HF_number,1:x_HF_number,len_index)*lenD(len_index);
                    dK_dcov{len_index} = dK_dlenD;
                end

                % eta D
                dK_dcov{vari_num+1} = KD;

                % len L
                for len_index = 1:vari_num
                    dK_dlenL = -KL.*sq_dis_v(:,:,len_index)*lenL(len_index);
                    dK_dcov{(vari_num+1)+len_index} = dK_dlenL;
                end

                % eta L
                dK_dcov{2*(vari_num+1)} = KL;

                % rho
                dK_drho = zeros(x_num);
                dK_drho(1:x_HF_number,1:x_HF_number) = ...
                    2*rho*rho*eta_exp_disL(1:x_HF_number,1:x_HF_number);
                dK_drho(1:x_HF_number,(x_HF_number+1):end) = ...
                    rho*eta_exp_disL(1:x_HF_number,(x_HF_number+1):end);
                dK_drho((x_HF_number+1):end,1:x_HF_number) = ...
                    dK_drho(1:x_HF_number,(x_HF_number+1):end)';
                dK_dcov{end} = dK_drho;
            end
        end

    end

    function [post,nlZ,dnlZ] = infEP(hyp,mean,cov,lik,x,y)
        % Expectation Propagation approximation to the posterior Gaussian Process.
        % The function takes a specified covariance function (see covFunctions.m) and
        % likelihood function (see likFunctions.m),and is designed to be used with
        % gp.m. See also infMethods.m. In the EP algorithm,the sites are
        % updated in random order,for better performance when cases are ordered
        % according to the targets.
        %
        % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-09-13.
        %
        % See also INFMETHODS.M.
        %
        persistent last_ttau last_tnu              % keep tilde parameters between calls
        tol = 1e-4; max_sweep = 10; min_sweep = 2;     % tolerance to stop EP iterations

        inf = 'infEP';
        if isnumeric(cov),K = cov;                    % use provided covariance matrix
        else K = feval(cov{:},hyp.cov,x); end       % evaluate the covariance matrix

        n = size(x,1);

        if isnumeric(mean),m = mean;                         % use provided mean vector
        else m = feval(mean{:},hyp.mean,x); end             % evaluate the mean vector

        % A note on naming: variables are given short but descriptive names in
        % accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
        % and s2 are mean and variance,nu and tau are natural parameters. A leading t
        % means tilde,a subscript _ni means "not i" (for cavity parameters),or _n
        % for a vector of cavity parameters. N(f|mu,Sigma) is the posterior.

        % marginal likelihood for ttau = tnu = zeros(n,1); equals n*log(2) for likCum*
        nlZ0 = -sum(feval(lik{:},hyp.lik,y,m,diag(K),inf));
        if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
            ttau = zeros(n,1); tnu  = zeros(n,1);        % init to zero if no better guess
            Sigma = K;                     % initialize Sigma and mu,the parameters of ..
            mu = m; nlZ = nlZ0;                  % .. the Gaussian posterior approximation
        else
            ttau = last_ttau; tnu  = last_tnu;   % try the tilde values from previous call
            [Sigma,mu,L,alpha,nlZ] = epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf);
            if nlZ > nlZ0                                           % if zero is better ..
                ttau = zeros(n,1); tnu  = zeros(n,1);       % .. then init with zero instead
                Sigma = K;                   % initialize Sigma and mu,the parameters of ..
                mu = m; nlZ = nlZ0;                % .. the Gaussian posterior approximation
            end
        end

        nlZ_old = Inf; sweep = 0;               % converged,max. sweeps or min. sweeps?
        while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep
            nlZ_old = nlZ; sweep = sweep+1;
            for i = randperm(n)       % iterate EP updates (in random order) over examples
                tau_ni = 1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
                nu_ni = mu(i)/Sigma(i,i)-tnu(i);                % .. params tau_ni and nu_ni

                % compute the desired derivatives of the indivdual log partition function
                [lZ,dlZ,d2lZ] = feval(lik{:},hyp.lik,y(i),nu_ni/tau_ni,1/tau_ni,inf);
                ttau_old = ttau(i); tnu_old = tnu(i);  % find the new tilde params,keep old
                ttau(i) = -d2lZ /(1+d2lZ/tau_ni);
                ttau(i) = max(ttau(i),0); % enforce positivity i.e. lower bound ttau by zero
                tnu(i)  = ( dlZ - nu_ni/tau_ni*d2lZ )/(1+d2lZ/tau_ni);

                dtt = ttau(i)-ttau_old; dtn = tnu(i)-tnu_old;      % rank-1 update Sigma ..
                si = Sigma(:,i); ci = dtt/(1+dtt*si(i));
                Sigma = Sigma - ci*si*si';                         % takes 70% of total time
                mu = mu - (ci*(mu(i)+si(i)*dtn)-dtn)*si;               % .. and recompute mu
            end
            % recompute since repeated rank-one updates can destroy numerical precision
            [Sigma,mu,L,alpha,nlZ] = epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf);
        end

        if sweep  ==  max_sweep && abs(nlZ-nlZ_old) > tol
            error('maximum number of sweeps exceeded in function infEP')
        end

        last_ttau = ttau; last_tnu = tnu;                       % remember for next call
        post.alpha = alpha; post.sW = sqrt(ttau); post.L = L;  % return posterior params

        if nargout>2                                           % do we want derivatives?
            dnlZ = hyp;                                   % allocate space for derivatives
            tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
            nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
            sW = sqrt(ttau);
            F = alpha*alpha'-repmat(sW,1,n).*(L\(L'\diag(sW)));   % covariance hypers
            [K,dK] = feval(cov{:},hyp.cov,x,[]);
            for i = 1:length(hyp.cov)
                dnlZ.cov(i) = -sum(sum(F.*dK{i}))/2;
            end
            for i = 1:numel(hyp.lik)                                   % likelihood hypers
                dlik = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf,i);
                dnlZ.lik(i) = -sum(dlik);
            end
            [junk,dlZ] = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf);% mean hyps
            for i = 1:numel(hyp.mean)
                dm = feval(mean{:},hyp.mean,x,i);
                dnlZ.mean(i) = -dlZ'*dm;
            end
        end
    end

    function [Sigma,mu,L,alpha,nlZ] = epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf)
        % function to compute the parameters of the Gaussian approximation,Sigma and
        % mu,and the negative log marginal likelihood,nlZ,from the current site
        % parameters,ttau and tnu. Also returns L (useful for predictions).
        %
        n = length(y);                                      % number of training cases
        sW = sqrt(ttau);                                        % compute Sigma and mu
        L = chol(eye(n)+sW*sW'.*K);                            % L'*L = B = eye(n)+sW*K*sW
        V = L'\(repmat(sW,1,n).*K);
        Sigma = K - V'*V;
        alpha = tnu-sW.*(L\(L'\(sW.*(K*tnu+m))));
        mu = K*alpha+m; v = diag(Sigma);

        tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
        nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
        lZ = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf);
        p = tnu-m.*ttau; q = nu_n-m.*tau_n;                        % auxiliary vectors
        nlZ = sum(log(diag(L))) - sum(lZ) - p'*Sigma*p/2 + (v'*p.^2)/2 ...
            - q'*((ttau./tau_n.*q-2*p).*v)/2 - sum(log(1+ttau./tau_n))/2;

    end

    function A = meanConst(hyp,x,i)
        % Constant mean function. The mean function is parameterized as:
        %
        % m(x) = c
        %
        % The hyperparameter is:
        %
        % hyp = [ c ]
        %
        % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2010-08-04.
        %
        % See also MEANFUNCTIONS.M.
        %
        if nargin<2,A = '1'; return; end             % report number of hyperparameters
        if numel(hyp)~=1,error('Exactly one hyperparameter needed.'),end
        c = hyp;
        if nargin == 2
            A = c*ones(size(x,1),1);                                       % evaluate mean
        else
            if i == 1
                A = ones(size(x,1),1);                                          % derivative
            else
                A = zeros(size(x,1),1);
            end
        end
    end

    function [varargout] = likErf(hyp,y,mu,s2,inf)
        % likErf - Error function or cumulative Gaussian likelihood function for binary
        % classification or probit regression. The expression for the likelihood is
        %   likErf(t) = (1+erf(t/sqrt(2)))/2 = normcdf(t).
        %
        % Several modes are provided,for computing likelihoods,derivatives and moments
        % respectively,see likFunctions.m for the details. In general,care is taken
        % to avoid numerical issues when the arguments are extreme.
        %
        % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2014-03-19.
        %
        % See also LIKFUNCTIONS.M.
        %
        if nargin<3,varargout = {'0'}; return; end   % report number of hyperparameters
        if nargin>1,y = sign(y); y(y == 0) = 1; else y = 1; end % allow only +/- 1 values
        if numel(y) == 0,y = 1; end

        if nargin<5                              % prediction mode if inf is not present
            y = y.*ones(size(mu));                                       % make y a vector
            s2zero = 1; if nargin>3&&numel(s2)>0&&norm(s2)>eps,s2zero = 0; end  % s2 == 0 ?
            if s2zero                                         % log probability evaluation
                lp = logphi(y.*mu);
            else                                                              % prediction
                lp = likErf(hyp,y,mu,s2,'infEP');
            end
            p = exp(lp); ymu = {}; ys2 = {};
            if nargout>1
                ymu = 2*p-1;                                                % first y moment
                if nargout>2
                    ys2 = 4*p.*(1-p);                                        % second y moment
                end
            end
            varargout = {lp,ymu,ys2};
        else                                                            % inference mode
            switch inf
                case 'infLaplace'
                    if nargin<6                                             % no derivative mode
                        f = mu; yf = y.*f;                            % product latents and labels
                        varargout = cell(nargout,1); [varargout{:}] = logphi(yf);   % query logphi
                        if nargout>1
                            varargout{2} = y.*varargout{2};
                            if nargout>3,varargout{4} = y.*varargout{4}; end
                        end
                    else                                                       % derivative mode
                        varargout = {[],[],[]};                         % derivative w.r.t. hypers
                    end

                case 'infEP'
                    if nargin<6                                             % no derivative mode
                        z = mu./sqrt(1+s2); dlZ = {}; d2lZ = {};
                        if numel(y) > 0,z = z.*y; end
                        if nargout <= 1,lZ = logphi(z);                         % log part function
                        else          [lZ,n_p] = logphi(z); end
                        if nargout > 1
                            if numel(y) == 0,y = 1; end
                            dlZ = y.*n_p./sqrt(1+s2);                      % 1st derivative wrt mean
                            if nargout>2,d2lZ = -n_p.*(z+n_p)./(1+s2); end         % 2nd derivative
                        end
                        varargout = {lZ,dlZ,d2lZ};
                    else                                                       % derivative mode
                        varargout = {[]};                                     % deriv. wrt hyp.lik
                    end
            end
        end
    end

    function [lp,dlp,d2lp,d3lp] = logphi(z)
        % Safe computation of logphi(z) = log(normcdf(z)) and its derivatives
        %                    dlogphi(z) = normpdf(x)/normcdf(x).
        % The function is based on index 5725 in Hart et al. and gsl_sf_log_erfc_e.
        %
        % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2013-11-13.
        %
        z = real(z);                                 % support for real arguments only
        lp = zeros(size(z));                                         % allocate memory
        id1 = z.*z<0.0492;                                 % first case: close to zero
        lp0 = -z(id1)/sqrt(2*pi);
        c = [ 0.00048204; -0.00142906; 0.0013200243174; 0.0009461589032;
            -0.0045563339802; 0.00556964649138; 0.00125993961762116;
            -0.01621575378835404; 0.02629651521057465; -0.001829764677455021;
            2*(1-pi/3); (4-pi)/3; 1; 1];
        f = 0; for i = 1:14,f = lp0.*(c(i)+f); end,lp(id1) = -2*f-log(2);
        id2 = z<-11.3137;                                    % second case: very small
        r = [ 1.2753666447299659525; 5.019049726784267463450;
            6.1602098531096305441; 7.409740605964741794425;
            2.9788656263939928886 ];
        q = [ 2.260528520767326969592;  9.3960340162350541504;
            12.048951927855129036034; 17.081440747466004316;
            9.608965327192787870698;  3.3690752069827527677 ];
        num = 0.5641895835477550741; for i = 1:5,num = -z(id2).*num/sqrt(2) + r(i); end
        den = 1.0;                   for i = 1:6,den = -z(id2).*den/sqrt(2) + q(i); end
        e = num./den; lp(id2) = log(e/2) - z(id2).^2/2;
        id3 = ~id2 & ~id1; lp(id3) = log(erfc(-z(id3)/sqrt(2))/2);  % third case: rest
        if nargout>1                                        % compute first derivative
            dlp = zeros(size(z));                                      % allocate memory
            dlp( id2) = abs(den./num) * sqrt(2/pi); % strictly positive first derivative
            dlp(~id2) = exp(-z(~id2).*z(~id2)/2-lp(~id2))/sqrt(2*pi); % safe computation
            if nargout>2                                     % compute second derivative
                d2lp = -dlp.*abs(z+dlp);             % strictly negative second derivative
                if nargout>3                                    % compute third derivative
                    d3lp = -d2lp.*abs(z+2*dlp)-dlp;     % strictly positive third derivative
                end
            end
        end
    end

end

%% surrogate model
function [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output] = getMFRBFModel...
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

function [predict_function_MFRBF, MFRBF_model] = interpMultiRadialBasisPreModel...
    (XHF, YHF, varargin)
% multi fildelity radial basis function interp pre model function
% XHF, YHF are x_HF_number x variable_number matrix
% XLF, YLF are x_LF_number x variable_number matrix
% aver_X, stdD_X is 1 x x_HF_number matrix
%
% input:
% XHF, YHF, basis_func_HF(can be []), XLF, YLF, basis_func_LF(can be [])
% XHF, YHF, basis_func_HF(can be []), LF_model
%
% output:
% predict_function, HK_model
%
% reference: [1] LIU Y, WANG S, ZHOU Q, et al. Modified Multifidelity
% Surrogate Model Based on Radial Basis Function with Adaptive Scale Factor
% [J]. Chinese Journal of Mechanical Engineering, 2022, 35(1): 77.
%
% Copyright 2023 Adel
%
[x_HF_number, variable_number] = size(XHF);
switch nargin
    case 4
        basis_func_HF = varargin{1};
        LF_model = varargin{2};

        % check whether LF model exist predict_function
        if ~isfield(LF_model, 'predict_function')
            error('interpHieraKrigingPreModel: low fidelity lack predict function');
        end
    case 6
        basis_func_HF = varargin{1};
        XLF = varargin{2};
        YLF = varargin{3};
        basis_func_LF = varargin{4};

        [x_LF_number, variable_number] = size(XLF);

        % first step
        % construct low fidelity model

        % normalize data
        aver_X = mean(XLF);
        stdD_X = std(XLF);
        aver_Y = mean(YLF);
        stdD_Y = std(YLF);
        index__ = find(stdD_X == 0);
        if  ~isempty(index__), stdD_X(index__) = 1; end
        index__ = find(stdD_Y == 0);
        if  ~isempty(index__), stdD_Y(index__) = 1; end

        XLF_nomlz = (XLF-aver_X)./stdD_X;
        YLF_nomlz = (YLF-aver_Y)./stdD_Y;

        if isempty(basis_func_LF)
            c = (prod(max(XLF_nomlz) - min(XLF_nomlz))/x_LF_number)^(1/variable_number);
            basis_func_LF = @(r) exp(-(r.^2)/c);
        end

        % initialization distance of XLF_nomlz
        XLF_dis = zeros(x_LF_number, x_LF_number);
        for variable_index = 1:variable_number
            XLF_dis = XLF_dis + ...
                (XLF_nomlz(:, variable_index) - XLF_nomlz(:, variable_index)').^2;
        end
        XLF_dis = sqrt(XLF_dis);

        [beta_LF, rdibas_matrix_LF] = interpRadialBasis...
            (XLF_dis, YLF_nomlz, basis_func_LF, x_LF_number);

        % initialization predict function
        predict_function_LF = @(X_predict) interpRadialBasisPredictor...
            (X_predict, XLF_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_LF_number, variable_number, beta_LF, basis_func_LF);

        LF_model.X = XLF;
        LF_model.Y = YLF;
        LF_model.radialbasis_matrix = rdibas_matrix_LF;
        LF_model.beta = beta_LF;

        LF_model.aver_X = aver_X;
        LF_model.stdD_X = stdD_X;
        LF_model.aver_Y = aver_Y;
        LF_model.stdD_Y = stdD_Y;
        LF_model.basis_function = basis_func_LF;

        LF_model.predict_function = predict_function_LF;
    otherwise
        error('interpHieraKrigingPreModel: error input');
end
MFRBF_model.LF_model = LF_model;
predict_function_LF = LF_model.predict_function;

% predict LF value at XHF point
YHF_pred = predict_function_LF(XHF);

% nomalizae
YHF_pred_nomlz = (YHF_pred - aver_Y)./stdD_Y;

% second step
% construct MFRBF model

% normalize data
aver_X = mean(XHF);
stdD_X = std(XHF);
aver_Y = mean(YHF);
stdD_Y = std(YHF);
index__ = find(stdD_X == 0);
if ~isempty(index__), stdD_X(index__) = 1;end
index__ = find(stdD_Y == 0);
if ~isempty(index__), stdD_Y(index__) = 1;end
XHF_nomlz = (XHF - aver_X)./stdD_X;
YHF_nomlz = (YHF - aver_Y)./stdD_Y;

if isempty(basis_func_HF)
    c = (prod(max(XHF_nomlz) - min(XHF_nomlz))/x_HF_number)^(1/variable_number);
    basis_func_HF = @(r) exp(-(r.^2)/c);
end

% initialization distance of XHF_nomlz
XHF_dis = zeros(x_HF_number, x_HF_number);
for variable_index = 1:variable_number
    XHF_dis = XHF_dis + ...
        (XHF_nomlz(:, variable_index) - XHF_nomlz(:, variable_index)').^2;
end
XHF_dis = sqrt(XHF_dis);

[beta_HF, rdibas_matrix_HF] = interpMultiRadialBasis...
    (XHF_dis, YHF_nomlz, basis_func_HF, x_HF_number, YHF_pred_nomlz);

% initialization predict function
predict_function_MFRBF = @(X_predict) interpMultiRadialBasisPredictor...
    (X_predict, XHF_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
    x_HF_number, variable_number, beta_HF, basis_func_HF, predict_function_LF);

MFRBF_model.X = XHF;
MFRBF_model.Y = YHF;
MFRBF_model.radialbasis_matrix = rdibas_matrix_HF;
MFRBF_model.beta = beta_HF;

MFRBF_model.aver_X = aver_X;
MFRBF_model.stdD_X = stdD_X;
MFRBF_model.aver_Y = aver_Y;
MFRBF_model.stdD_Y = stdD_Y;
MFRBF_model.basis_function = basis_func_HF;

MFRBF_model.predict_function = predict_function_MFRBF;

% abbreviation:
% num: number, pred: predict, vari: variable
    function [beta, rdibas_matrix] = interpRadialBasis...
            (X_dis, Y, basis_function, x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix = rdibas_matrix + eye(x_number)*1e-6;

        % solve beta
        beta = rdibas_matrix\Y;
    end

    function [Y_pred] = interpRadialBasisPredictor...
            (X_pred, X_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_num, vari_num, beta, basis_function)
        % radial basis function interpolation predict function
        %
        [x_pred_num, ~] = size(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred - aver_X)./stdD_X;

        % calculate distance
        X_dis_pred = zeros(x_pred_num, x_num);
        for vari_index = 1:vari_num
            X_dis_pred = X_dis_pred + ...
                (X_pred_nomlz(:, vari_index) - X_nomlz(:, vari_index)').^2;
        end
        X_dis_pred = sqrt(X_dis_pred);

        % predict variance
        Y_pred = basis_function(X_dis_pred)*beta;

        % normalize data
        Y_pred = Y_pred*stdD_Y + aver_Y;
    end

    function [beta, rdibas_matrix] = interpMultiRadialBasis...
            (X_dis, Y, basis_function, x_number, YHF_pred)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix = rdibas_matrix + eye(x_number)*1e-6;

        % add low fildelity value
        H = [rdibas_matrix.*YHF_pred, rdibas_matrix];

        % solve beta
        beta = H'*((H*H')\Y);
    end

    function [Y_pred] = interpMultiRadialBasisPredictor...
            (X_pred, X_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_num, vari_num, beta, basis_function, predict_function_LF)
        % radial basis function interpolation predict function
        %
        [x_pred_num, ~] = size(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred - aver_X)./stdD_X;

        % calculate distance
        X_dis_pred = zeros(x_pred_num, x_num);
        for vari_index = 1:vari_num
            X_dis_pred = X_dis_pred + ...
                (X_pred_nomlz(:, vari_index) - X_nomlz(:, vari_index)').^2;
        end
        X_dis_pred = sqrt(X_dis_pred);

        % predict low fildelity value
        Y_pred_LF = predict_function_LF(X_pred);

        % nomalizae
        Y_pred_LF_nomlz = (Y_pred_LF - aver_Y)./stdD_Y;

        % combine two matrix
        rdibas_matrix_pred = basis_function(X_dis_pred);
        H = [rdibas_matrix_pred.*Y_pred_LF_nomlz, rdibas_matrix_pred];

        % predict variance
        Y_pred = H*beta;

        % normalize data
        Y_pred = Y_pred*stdD_Y + aver_Y;
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
