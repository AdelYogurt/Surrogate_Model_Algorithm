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
% model_function = modelFunction(x,@(x) benchmark.singleGPObject(x),[]);

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
% model_function = modelFunction(x,@(x) benchmark.singlePKObject(x),[]);

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
model_function = @(x) modelFunction(x,@(x) benchmark.singleEP20Object(x),[]);

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
% nonlcon_function = @(x) benchmark.singlePVD4Nonlcon(x);
% model_function = @(x) modelFunction(x,@(x) benchmark.singlePVD4Object(x),@(x) violationFunction(x,A,B,Aeq,Beq,@(x) benchmark.singlePVD4Nonlcon(x)));
% cheapcon_function = [];

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
% nonlcon_function = [];
% nonlcon_function_LF = [];
% model_function = @(x) modelFunction(x,@(x) benchmark.singleG01Object(x),@(x) violationFunction(x,A,B,Aeq,Beq,[]));
% cheapcon_function = [];

% x_initial = rand(1,variable_number).*(up_bou-low_bou)+low_bou;
% [x_best,fval_best,~,output] = fmincon(object_function,x_initial,A,B,Aeq,Beq,low_bou,up_bou,[],optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',10000,'Display','iter-detailed'))

%% single run

delete([data_library_name,'.txt']);
delete('result_total.txt');

[x_best,fval_best,NFE,output] = optimalSurrogateSADE...
    (model_function,variable_number,low_bou,up_bou,...
    cheapcon_function,200,300)
result_x_best = output.result_x_best;
result_fval_best = output.result_fval_best;

figure(1);
plot(result_fval_best);

figure(2);
[x_list,fval_list,con_list,coneq_list] = dataLibraryRead...
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
%     [x_best,fval_best,NFE,output] = optimalSurrogateSADE...
%         (model_function,variable_number,low_bou,up_bou,...
%         cheapcon_function,max_NFE,300);
%
%     result_fval(repeat_index) = fval_best;
% end
%
% fprintf('Fval     : lowest = %4.4f,mean = %4.4f,worst = %4.4f,std = %4.4f \n',min(result_fval),mean(result_fval),max(result_fval),std(result_fval));
% object_function_name = char(object_function);
% save([object_function_name(15:end-3),'_',num2str(max_NFE),'_SADE','.mat']);

%% main
function [x_best,fval_best,NFE,output] = optimalSurrogateSADE...
    (model_function,variable_number,low_bou,up_bou,...
    cheapcon_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance,x_initial_list)
% KRG-CDE optimization algorithm
%
% referance: [1] 叶年辉,龙腾,武宇飞,et al.
% 基于Kriging代理模型的约束差分进化算法 [J]. 航空学报,2021,42(6): 13.
%
% Copyright 2022 Adel
%
if nargin < 10
    x_initial_list = [];
    if nargin < 9 || isempty(nonlcon_torlance)
        nonlcon_torlance = 1e-3;
        if nargin < 8 || isempty(torlance)
            torlance = 1e-3;
            if nargin < 7
                iteration_max = [];
                if nargin < 6
                    NFE_max = [];
                end
            end
        end
    end
end

if nargin < 5
    cheapcon_function = [];
end

DRAW_FIGURE_FLAG = 0; % whether draw data
INFORMATION_FLAG = 1; % whether print data
CONVERGENCE_JUDGMENT_FLAG = 0; % whether judgment convergence

% hyper parameter
population_number = min(100,10*variable_number);
RBF_number = max(100,(variable_number+1)*(variable_number+2)/2);
scaling_factor = 0.8; % F
cross_rate = 0.8;

% max fval when normalize fval,con,coneq
nomlz_fval = 100;

protect_range = 1e-5;

data_library_name = 'optimal_data_library';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done = 0;NFE = 0;iteration = 0;

% step 2
% generate initial sample x_list
if isempty(x_initial_list)
    %     [~,x_updata_list,~] = getLatinHypercube...
    %         (population_number,variable_number,[],low_bou,up_bou,cheapcon_function);
    x_updata_list = lhsdesign(population_number,variable_number).*(up_bou-low_bou)+low_bou;
else
    x_updata_list = x_initial_list;
end

% detech expensive constraints and initializa data library
if ~isempty(x_updata_list)
    [x_list,fval_list,con_list,coneq_list] = dataLibraryWrite...
        (data_library_name,model_function,x_updata_list(1,:));NFE = NFE+1;
    x_updata_list = x_updata_list(2:end,:);
else
    [x_list,fval_list,con_list,coneq_list] = dataLibraryRead(data_library_name,low_bou,up_bou);
end
vio_list = calViolation(con_list,coneq_list,nonlcon_torlance);
if ~isempty(vio_list)
    expensive_nonlcon_flag = 1;
else
    expensive_nonlcon_flag = 0;
end
data_library = struct('x_list',x_list,'fval_list',fval_list,'con_list',con_list,'coneq_list',coneq_list,'vio_list',vio_list);

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

% updata data library by x_list
[x_updata_list,fval_updata_list,con_updata_list,coneq_updata_list] = dataLibraryWrite...
    (data_library_name,model_function,x_updata_list);NFE = NFE+size(x_updata_list,1);
vio_updata_list = calViolation(con_updata_list,coneq_updata_list,nonlcon_torlance);
[data_library,x_list,fval_list,con_list,coneq_list,vio_list] = dataLibraryUpdata...
    (data_library,x_updata_list,fval_updata_list,con_updata_list,coneq_updata_list,vio_updata_list);

iteration = iteration+1;
kriging_model_fval = [];
kriging_model_con = [];
kriging_model_coneq = [];
next_search_flag = 'G'; % 'G' is global search,'l' is local search
while ~done
    search_flag = next_search_flag;
    % nomalization con with average
    fval_max = max(abs(fval_list),[],1);
    fval_nomlz_list = fval_list;
    if ~isempty(con_list)
        con_max_list = max(abs(con_list),[],1);
        con_nomlz_list = con_list;
    else
        con_nomlz_list = [];
    end
    if ~isempty(coneq_list)
        coneq_max_list = max(abs(coneq_list),[],1);
        coneq_nomlz_list = coneq_list;
    else
        coneq_nomlz_list = [];
    end
    vio_nomlz_list = calViolation(con_nomlz_list,coneq_nomlz_list,nonlcon_torlance);

    % find fesiable data in current data library
    if expensive_nonlcon_flag
        feasi_boolean_list = vio_list <= 0;
    end

    if search_flag == 'G'
        % global search
        [x_global_infill,...
            kriging_model_fval,kriging_model_con,kriging_model_coneq] = searchGlobal...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
            population_number,scaling_factor,cross_rate,...
            kriging_model_fval,kriging_model_con,kriging_model_coneq,...
            expensive_nonlcon_flag);

        [x_global_infill,fval_global_infill,con_global_infill,coneq_global_infill,NFE_p,repeat_index] = dataLibraryWriteProtect...
            (data_library_name,model_function,x_global_infill,...
            x_list,low_bou,up_bou,protect_range);NFE = NFE+NFE_p;

        % infill point violation and updata to library
        vio_global_infill = calViolation(con_global_infill,coneq_global_infill,nonlcon_torlance);
        [data_library,x_list,fval_list,con_list,coneq_list,vio_list] = dataLibraryUpdata...
            (data_library,x_global_infill,fval_global_infill,con_global_infill,coneq_global_infill,vio_global_infill);

        next_search_flag = 'l';
        if isempty(x_global_infill)
            % process error
            x_global_infill = x_list(repeat_index,:);
            fval_global_infill = fval_list(repeat_index,:);
            if ~isempty(con_list)
                con_global_infill = con_list(repeat_index,:);
            end
            if ~isempty(coneq_list)
                coneq_global_infill = coneq_list(repeat_index,:);
            end
            if ~isempty(vio_list)
                vio_global_infill = vio_list(repeat_index,:);
            end

        else
            if expensive_nonlcon_flag
                min_vio = min(vio_list);
                min_fval = min(fval_list([feasi_boolean_list;true(0)]),[],1);

                % if all point is infeasible,violation of point infilled is
                % less than min violation of all point means improve.if
                % feasible point exist,fval of point infilled is less than min
                % fval means improve
                if vio_global_infill < min_vio
                    if ~isempty(min_fval)
                        if fval_global_infill < min_fval
                            % improve, continue global search
                            next_search_flag = 'G';
                        end
                    else
                        next_search_flag = 'G';
                    end
                end
            else
                min_fval = min(fval_list(1:end-1));

                % imporve, continue global search
                if fval_global_infill < min_fval
                    next_search_flag = 'G';
                end
            end

            if DRAW_FIGURE_FLAG && variable_number < 3
                interpVisualize(kriging_model_fval,low_bou,up_bou);
                line(x_global_infill(1),x_global_infill(2),fval_global_infill./fval_max*nomlz_fval,'Marker','o','color','r');
            end

        end

    else
        % local search
        [x_local_infill,...
            RBF_model_fval,RBF_model_con,RBF_model_coneq] = searchLocal...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
            population_number,RBF_number,...
            expensive_nonlcon_flag);

        [x_local_infill,fval_local_infill,con_local_infill,coneq_local_infill,NFE_p,repeat_index] = dataLibraryWriteProtect...
            (data_library_name,model_function,x_local_infill,...
            x_list,low_bou,up_bou,protect_range);NFE = NFE+NFE_p;

        % infill point violation and updata to library
        vio_local_infill = calViolation(con_local_infill,coneq_local_infill,nonlcon_torlance);
        [data_library,x_list,fval_list,con_list,coneq_list,vio_list] = dataLibraryUpdata...
            (data_library,x_local_infill,fval_local_infill,con_local_infill,coneq_local_infill,vio_local_infill);

        next_search_flag = 'G';
        if isempty(x_local_infill)
            % process error
            x_local_infill = x_list(repeat_index,:);
            fval_local_infill = fval_list(repeat_index,:);
            if ~isempty(con_list)
                con_local_infill = con_list(repeat_index,:);
            end
            if ~isempty(coneq_list)
                coneq_local_infill = coneq_list(repeat_index,:);
            end
            if ~isempty(vio_list)
                vio_local_infill = vio_list(repeat_index,:);
            end
        else
            if expensive_nonlcon_flag
                min_vio = min(vio_list);
                min_fval = min(fval_list([feasi_boolean_list;true(0)]),[],1);

                % if all point is infeasible,violation of point infilled is
                % less than min violation of all point means improve.if
                % feasible point exist,fval of point infilled is less than min
                % fval means improve
                if vio_local_infill < min_vio
                    if ~isempty(min_fval)
                        if fval_local_infill < min_fval
                            % improve, continue local search
                            next_search_flag = 'l';
                        end
                    else
                        % improve, continue local search
                        next_search_flag = 'l';
                    end
                end
            else
                min_fval = min(fval_list(1:end-1));

                % fval of point infilled is less than min fval means improve
                if fval_local_infill < min_fval
                    % imporve, continue local search
                    next_search_flag = 'l';
                end
            end

            if DRAW_FIGURE_FLAG && variable_number < 3
                interpVisualize(RBF_model_fval,low_bou,up_bou);
                line(x_local_infill(1),x_local_infill(2),fval_local_infill./fval_max*nomlz_fval,'Marker','o','color','r');
            end
        end
    
    end

    % find best result to record
    [x_best,fval_best,con_best,coneq_best] = findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);
    vio_best = calViolation(con_best,coneq_best,nonlcon_torlance);

    if INFORMATION_FLAG
        fprintf('fval:    %f    violation:    %f    NFE:    %-3d\n',fval_best,vio_best,NFE);
        %         fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
        %         if search_flag == 0
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
end

result_x_best = result_x_best(1:iteration-1,:);
result_fval_best = result_fval_best(1:iteration-1);

output.result_x_best = result_x_best;
output.result_fval_best = result_fval_best;
output.data_library = data_library;

    function [x_updata_list,fval_updata_list,con_updata_list,coneq_updata_list,NFE_updata,repeat_index] = dataLibraryWriteProtect...
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
                [x_updata__,fval_updata__,con_updata__,coneq_updata__] = dataLibraryWrite...
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
function vio_list = calViolation(con_list,coneq_list,nonlcon_torlance)
% calculate violation of data
%
if isempty(con_list) && isempty(coneq_list)
    vio_list = [];
else
    vio_list = zeros(max(size(con_list,1),size(coneq_list,1)),1);
    if ~isempty(con_list)
        vio_list = vio_list+sum(max(con_list-nonlcon_torlance,0),2);
    end
    if ~isempty(coneq_list)
        vio_list = vio_list+sum((abs(coneq_list)-nonlcon_torlance),2);
    end

end
end

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

function [x_list,fval_list,con_list,coneq_list,vio_list] = rankData...
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

function [x_global_infill,...
    kriging_model_fval,kriging_model_con,kriging_model_coneq] = searchGlobal...
    (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
    variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
    population_number,scaling_factor,cross_rate,...
    kriging_model_fval,kriging_model_con,kriging_model_coneq,...
    expensive_nonlcon_flag)
% find global infill point function
%

% step 5
% rank x_list data
[x_rank_list,~,~,~] = rankData...
    (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
    cheapcon_function,nonlcon_torlance);

% step 6
% only the first population_number will be use
x_best_popu_list = x_rank_list(1:population_number,:);

% differ evolution mutations
X_new_R1 = differEvolutionRand...
    (low_bou,up_bou,x_best_popu_list,scaling_factor,population_number,1);
X_new_R2 = differEvolutionRand...
    (low_bou,up_bou,x_best_popu_list,scaling_factor,population_number,2);
X_new_CR = differEvolutionCurrentRand...
    (low_bou,up_bou,x_best_popu_list,scaling_factor);
X_new_CB = differEvolutionCurrentBest...
    (low_bou,up_bou,x_best_popu_list,scaling_factor,1);

% differ evolution crossover
X_new_R1 = differEvolutionCrossover...
    (low_bou,up_bou,x_best_popu_list,X_new_R1,cross_rate);
X_new_R2 = differEvolutionCrossover...
    (low_bou,up_bou,x_best_popu_list,X_new_R2,cross_rate);
X_new_CR = differEvolutionCrossover...
    (low_bou,up_bou,x_best_popu_list,X_new_CR,cross_rate);
X_new_CB = differEvolutionCrossover...
    (low_bou,up_bou,x_best_popu_list,X_new_CB,cross_rate);

% find global infill point base kriging model from offspring X
x_DE_list = [X_new_R1;X_new_R2;X_new_CR;X_new_CB];

% step 4
% updata kriging model and function
[kriging_model_fval,kriging_model_con,kriging_model_coneq,output_kriging] = getKrigingModel...
    (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
    kriging_model_fval,kriging_model_con,kriging_model_coneq);
object_function_surrogate = output_kriging.object_function_surrogate;
nonlcon_function_surrogate = output_kriging.nonlcon_function_surrogate;

%         % modify
%         srgtOPT     =  srgtsKRGSetOptions(x_list , fval_nomlz_list);
%         srgt_KRG    = srgtsKRGFit(srgtOPT);
%         [fval_pred_DE_list, fval_var_DE_list] = srgtsKRGPredictor(x_DE_list, srgt_KRG);

% evaluate each x_offspring fval and constraints
[fval_pred_DE_list,fval_var_DE_list] = object_function_surrogate(x_DE_list);
if expensive_nonlcon_flag

    if ~isempty(nonlcon_function_surrogate)
        [con_pred_DE_list,con_var_DE_list,coneq_pred_DE_list,coneq_var_DE_list] = nonlcon_function_surrogate(x_DE_list);
    end

    %             % modify
    %             srgtOPT     =  srgtsKRGSetOptions(x_list , con_nomlz_list);
    %             srgt_KRG    = srgtsKRGFit(srgtOPT);
    %             [con_pred_DE_list, con_var_DE_list] = srgtsKRGPredictor(x_DE_list, srgt_KRG);

    vio_DE_list = zeros(4*population_number,1);
    if ~isempty(con_nomlz_list)
        vio_DE_list = vio_DE_list+sum(max(con_pred_DE_list-nonlcon_torlance,0),2);
    end
    if ~isempty(coneq_nomlz_list)
        vio_DE_list = vio_DE_list+sum((abs(con_pred_DE_list)-nonlcon_torlance),2);
    end
    feasi_boolean_DE_list = vio_DE_list <= nonlcon_torlance;
else
    feasi_boolean_DE_list = true(ones(1,4*population_number));
end

% if have feasiable_index_list,only use feasiable to choose
if all(~feasi_boolean_DE_list)
    % base on constaints improve select global infill
    % lack process of equal constraints
    con_nomlz_base = max(min(con_nomlz_list,[],1),0);
    con_impove_probability_list = sum(...
        normcdf((con_nomlz_base-con_pred_DE_list)./sqrt(con_var_DE_list)),2);
    [~,con_best_index] = max(con_impove_probability_list);
    con_best_index = con_best_index(1);
    x_global_infill = x_DE_list(con_best_index,:);
else
    % base on fitness DE point to select global infill
    if expensive_nonlcon_flag
        x_DE_list = x_DE_list(feasi_boolean_DE_list,:);
        fval_pred_DE_list = fval_pred_DE_list(feasi_boolean_DE_list);
        fval_var_DE_list = fval_var_DE_list(feasi_boolean_DE_list);
    end

    fval_DE_min = min(fval_pred_DE_list,[],1);
    fval_DE_max = max(fval_pred_DE_list,[],1);
    fval_var_DE_min = min(fval_var_DE_list,[],1);
    fval_var_DE_max = max(fval_var_DE_list,[],1);
    % modify
    %             DE_fitness_list = -fval_DE_list+fval_var_DE_list;
    DE_fitness_list = -(fval_pred_DE_list-fval_DE_min)/(fval_DE_max-fval_DE_min)+...
        (fval_var_DE_list-fval_var_DE_min)/(fval_var_DE_max-fval_var_DE_min);
    [~,fitness_best_index] = max(DE_fitness_list);
    fitness_best_index = fitness_best_index(1);
    x_global_infill = x_DE_list(fitness_best_index,:);
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

end

function [x_local_infill,...
    RBF_model_fval,RBF_model_con,RBF_model_coneq] = searchLocal...
    (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
    variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
    population_number,RBF_number,...
    expensive_nonlcon_flag)
% find local infill point function
%

[x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,vio_nomlz_list] = rankData...
    (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
    cheapcon_function,nonlcon_torlance);

% step 8
% rand select initial local point from x_list
x_index = randi(population_number);
%         x_index = find(vio_nomlz_list == 0);
%         x_index = x_index(end)+1;
x_initial = x_list(x_index,:);

%         % rank x_list data and select best point as initial local point
%         [x_list_rank,~,~,~] = rankData...
%             (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
%             cheapcon_function,nonlcon_torlance);
%         x_initial = x_list_rank(1,:);

% select nearest point to construct RBF
RBF_number = min(RBF_number,size(x_list,1));
distance = sum(((x_initial-x_list)./(up_bou-low_bou)).^2,2);
[~,index_list] = sort(distance);
index_list = index_list(1:RBF_number);
x_RBF_list = x_list(index_list,:);
fval_RBF_nomlz_list = fval_nomlz_list(index_list,:);
if ~isempty(con_nomlz_list)
    con_RBF_nomlz_list = con_nomlz_list(index_list,:);
else
    con_RBF_nomlz_list = [];
end
if ~isempty(coneq_nomlz_list)
    coneq_RBF_nomlz_list = coneq_nomlz_list(index_list,:);
else
    coneq_RBF_nomlz_list = [];
end

% modify
% get RBF model and function
[RBF_model_fval,RBF_model_con,RBF_model_coneq,output_RBF] = getRadialBasisModel...
    (x_RBF_list,fval_RBF_nomlz_list,con_RBF_nomlz_list,coneq_RBF_nomlz_list);
object_function_surrogate = output_RBF.object_function_surrogate;
nonlcon_function_surrogate = output_RBF.nonlcon_function_surrogate;
low_bou_local = min(x_RBF_list,[],1);
up_bou_local = max(x_RBF_list,[],1);

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
    low_bou_local,up_bou_local,constraint_function,fmincon_options);

    function Obj = SrgtObj(SRGTRBF_structObj,x)
        [Obj, ~] = RBFPredict(x, SRGTRBF_structObj);
        % Obj = Obj + pred;
    end

    function [g, h] = SrgtCon(SRGTRBF_structCong,x)
        h = [];
        g = RBFPredict(x, SRGTRBF_structCong);
    end

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

%% surrogate model
function [kriging_model_fval, kriging_model_con, kriging_model_coneq, output] = getKrigingModel...
    (x_list, fval_list, con_list, coneq_list, ...
    kriging_model_fval, kriging_model_con, kriging_model_coneq)
% base on library_data to create kriging model and function
% if input model, function will updata model
% object_function is multi fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
if nargin < 5
    kriging_model_fval = [];
end

if size(x_list, 1) ~= size(fval_list, 1)
    error('getKrigingModel: x_list size no equal fval_list size')
end

if isempty(kriging_model_fval)
    [predict_function_fval, kriging_model_fval] = interpKrigingPreModel...
        (x_list, fval_list);
else
    [predict_function_fval, kriging_model_fval] = interpKrigingPreModel...
        (x_list, fval_list, kriging_model_fval.hyp);
end

if ~isempty(con_list)
    predict_function_con = cell(size(con_list, 2), 1);
    if size(x_list, 1) ~= size(con_list, 1)
        error('getKrigingModel: x_list size no equal con_list size')
    end
    if isempty(kriging_model_con)
        kriging_model_con = struct('X', [], 'Y', [], ...
            'fval_regression', [], 'covariance', [], 'inv_covariance', [], ...
            'hyp', [], 'beta', [], 'gama', [], 'sigma_sq', [], ...
            'aver_X', [], 'stdD_X', [], 'aver_Y', [], 'stdD_Y', [], ...
            'predict_function', []);
        kriging_model_con = repmat(kriging_model_con, 1, [size(con_list, 2)]);
        for con_index = 1:size(con_list, 2)
            [predict_function_con{con_index}, kriging_model_con(con_index)] = interpKrigingPreModel...
                (x_list, con_list(:, con_index));
        end
    else
        for con_index = 1:size(con_list, 2)
            [predict_function_con{con_index}, kriging_model_con(con_index)] = interpKrigingPreModel...
                (x_list, con_list(:, con_index), kriging_model_con(con_index).hyp);
        end
    end
else
    predict_function_con = [];
    kriging_model_con = [];
end

if ~isempty(coneq_list)
    predict_function_coneq = cell(size(coneq_list, 2), 1);
    if size(x_list, 1) ~= size(coneq_list, 1)
        error('getKrigingModel: x_list size no equal coneq_list size')
    end
    if isempty(kriging_model_coneq)
        kriging_model_coneq = struct('X', [], 'Y', [], ...
            'fval_regression', [], 'covariance', [], 'inv_covariance', [], ...
            'hyp', [], 'beta', [], 'gama', [], 'sigma_sq', [], ...
            'aver_X', [], 'stdD_X', [], 'aver_Y', [], 'stdD_Y', [], ...
            'predict_function', []);
        kriging_model_coneq = repmat(kriging_model_coneq, 1, [size(coneq_list, 2)]);
        for coneq_index = 1:size(coneq_list, 2)
            [predict_function_coneq{coneq_index}, kriging_model_coneq(coneq_index)] = interpKrigingPreModel...
                (x_list, coneq_list(:, coneq_index));
        end
    else
        for coneq_index = 1:size(coneq_list, 2)
            [predict_function_coneq{coneq_index}, kriging_model_coneq(coneq_index)] = interpKrigingPreModel...
                (x_list, coneq_list(:, coneq_index), kriging_model_coneq(coneq_index).hyp);
        end
    end
else
    predict_function_coneq = [];
    kriging_model_coneq = [];
end

object_function_surrogate = @(X_predict) objectFunctionSurrogate(X_predict, predict_function_fval);
if isempty(predict_function_con) && isempty(predict_function_coneq)
    nonlcon_function_surrogate = [];
else
    nonlcon_function_surrogate = @(X_predict) nonlconFunctionSurrogate(X_predict, predict_function_con, predict_function_coneq);
end

output.object_function_surrogate = object_function_surrogate;
output.nonlcon_function_surrogate = nonlcon_function_surrogate;

    function [fval, fval_var] = objectFunctionSurrogate...
            (X_predict, predict_function_fval)
        % connect all predict favl
        %
        [fval, fval_var] = predict_function_fval(X_predict);
    end

    function [con, con_var, coneq, coneq_var] = nonlconFunctionSurrogate...
            (X_predict, predict_function_con, predict_function_coneq)
        % connect all predict con and coneq
        %
        if isempty(predict_function_con)
            con = [];
            con_var = [];
        else
            con = zeros(size(X_predict, 1), length(predict_function_con));
            con_var = zeros(size(X_predict, 1), length(predict_function_con));
            for con_index__ = 1:length(predict_function_con)
                [con(:, con_index__), con_var(:, con_index__)] = ....
                    predict_function_con{con_index__}(X_predict);
            end
        end
        if isempty(predict_function_coneq)
            coneq = [];
            coneq_var = [];
        else
            coneq = zeros(size(X_predict, 1), length(predict_function_coneq));
            coneq_var = zeros(size(X_predict, 1), length(predict_function_coneq));
            for coneq_index__ = 1:length(predict_function_coneq)
                [coneq(:, coneq_index__), coneq_var(:, coneq_index__)] = ...
                    predict_function_coneq{coneq_index__}(X_predict);
            end
        end
    end
end

function [predict_function, kriging_model] = interpKrigingPreModel...
    (X, Y, hyp)
% nomalization method is grassian
% add multi x_predict input support
% prepare model, optimal theta and calculation parameter
% X, Y are x_number x variable_number matrix
% aver_X, stdD_X is 1 x x_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
% theta = exp(hyp)
%
% input initial data X, Y, which are real data
%
% output is a kriging model, include predict_function...
% X, Y, base_function_list
%
% Copyright 2023.2 Adel
%
[x_number, variable_number] = size(X);
if nargin < 3
    hyp = 0;
end

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__ = find(stdD_X == 0);
if  ~isempty(index__), stdD_X(index__) = 1; end
index__ = find(stdD_Y == 0);
if  ~isempty(index__), stdD_Y(index__) = 1; end
X_nomlz = (X-aver_X)./stdD_X;
Y_nomlz = (Y-aver_Y)./stdD_Y;

% initial X_dis_sq
X_dis_sq = zeros(x_number, x_number, variable_number);
for variable_index = 1:variable_number
    X_dis_sq(:, :, variable_index) = ...
        (X_nomlz(:, variable_index)-X_nomlz(:, variable_index)').^2;
end

% regression function define
% notice reg_function process no normalization data
% reg_function = @(X) regZero(X);
reg_function = @(X) regLinear(X);

% calculate reg
fval_reg_nomlz = (reg_function(X)-aver_Y)./stdD_Y;

% optimal to get hyperparameter
fmincon_option = optimoptions('fmincon', 'Display', 'none', ...
    'OptimalityTolerance', 1e-2, ...
    'FiniteDifferenceStepSize', 1e-5, ...,
    'MaxIterations', 10, 'SpecifyObjectiveGradient', false);
low_bou_hyp = -3;
up_bou_hyp = 3;
object_function_hyp = @(hyp) objectNLLKriging...
    (X_dis_sq, Y_nomlz, x_number, variable_number, hyp, fval_reg_nomlz);

% [fval, gradient] = object_function_hyp(hyp)
% [~, gradient_differ] = differ(object_function_hyp, hyp)

% drawFunction(object_function_hyp, low_bou_hyp, up_bou_hyp);

hyp = fmincon...
    (object_function_hyp, hyp, [], [], [], [], low_bou_hyp, up_bou_hyp, [], fmincon_option);

% get parameter
[covariance, inv_covariance, beta, sigma_sq] = interpKriging...
    (X_dis_sq, Y_nomlz, x_number, variable_number, exp(hyp), fval_reg_nomlz);
gama = inv_covariance*(Y_nomlz-fval_reg_nomlz*beta);
FTRF = fval_reg_nomlz'*inv_covariance*fval_reg_nomlz;

% initialization predict function
predict_function = @(X_predict) interpKrigingPredictor...
    (X_predict, X_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
    x_number, variable_number, exp(hyp), beta, gama, sigma_sq, ...
    inv_covariance, fval_reg_nomlz, FTRF, reg_function);

kriging_model.X = X;
kriging_model.Y = Y;
kriging_model.fval_regression = fval_reg_nomlz;
kriging_model.covariance = covariance;
kriging_model.inv_covariance = inv_covariance;

kriging_model.hyp = hyp;
kriging_model.beta = beta;
kriging_model.gama = gama;
kriging_model.sigma_sq = sigma_sq;
kriging_model.aver_X = aver_X;
kriging_model.stdD_X = stdD_X;
kriging_model.aver_Y = aver_Y;
kriging_model.stdD_Y = stdD_Y;

kriging_model.predict_function = predict_function;

% abbreviation:
% num: number, pred: predict, vari: variable, hyp: hyper parameter
% NLL: negative log likelihood
    function [fval, gradient] = objectNLLKriging...
            (X_dis_sq, Y, x_num, vari_num, hyp, F_reg)
        % function to minimize sigma_sq
        %
        theta = exp(hyp);
        [cov, inv_cov, ~, sigma2, inv_FTRF, Y_Fmiu] = interpKriging...
            (X_dis_sq, Y, x_num, vari_num, theta, F_reg);

        % calculation negative log likelihood
        L = chol(cov)';
        fval = x_num/2*log(sigma2)+sum(log(diag(L)));

        % calculate gradient
        if nargout > 1
            % gradient
            dcov_dtheta = zeros(x_num, x_num);
            for vari_index = 1:vari_num
                dcov_dtheta = dcov_dtheta + X_dis_sq(:, :, vari_index);
            end
            dcov_dtheta = -dcov_dtheta.*cov*theta/vari_num;

            dinv_cov_dtheta = ...
                -inv_cov*dcov_dtheta*inv_cov;

            dinv_FTRF_dtheta = -inv_FTRF*...
                (F_reg'*dinv_cov_dtheta*F_reg)*...
                inv_FTRF;

            dmiu_dtheta = dinv_FTRF_dtheta*(F_reg'*inv_cov*Y)+...
                inv_FTRF*(F_reg'*dinv_cov_dtheta*Y);

            dY_Fmiu_dtheta = -F_reg*dmiu_dtheta;

            dsigma2_dtheta = (dY_Fmiu_dtheta'*inv_cov*Y_Fmiu+...
                Y_Fmiu'*dinv_cov_dtheta*Y_Fmiu+...
                Y_Fmiu'*inv_cov*dY_Fmiu_dtheta)/x_num;

            dlnsigma2_dtheta = 1/sigma2*dsigma2_dtheta;

            dlndetR = trace(inv_cov*dcov_dtheta);

            gradient = x_num/2*dlnsigma2_dtheta+0.5*dlndetR;

        end
    end

    function [cov, inv_cov, beta, sigma_sq, inv_FTRF, Y_Fmiu] = interpKriging...
            (X_dis_sq, Y, x_num, vari_num, theta, F_reg)
        % kriging interpolation kernel function
        % Y(x) = beta+Z(x)
        %
        cov = zeros(x_num, x_num);
        for vari_index = 1:vari_num
            cov = cov+X_dis_sq(:, :, vari_index)*theta;
        end
        cov = exp(-cov/vari_num)+eye(x_num)*1e-3;

        % coefficient calculation
        inv_cov = cov\eye(x_num);
        inv_FTRF = (F_reg'*inv_cov*F_reg)\eye(size(F_reg, 2));

        % basical bias
        beta = inv_FTRF*(F_reg'*inv_cov*Y);
        Y_Fmiu = Y-F_reg*beta;
        sigma_sq = (Y_Fmiu'*inv_cov*Y_Fmiu)/x_num;

    end

    function [Y_pred, Var_pred] = interpKrigingPredictor...
            (X_pred, X_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_num, vari_num, theta, beta, gama, sigma_sq, ...
            inv_cov, fval_reg_nomlz, FTRF, reg_function)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num, ~] = size(X_pred);
        fval_reg_pred = reg_function(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred-aver_X)./stdD_X;
        fval_reg_pred_nomlz = (fval_reg_pred-aver_Y)./stdD_Y;

        % predict covariance
        predict_cov = zeros(x_num, x_pred_num);
        for vari_index = 1:vari_num
            predict_cov = predict_cov+...
                (X_nomlz(:, vari_index)-X_pred_nomlz(:, vari_index)').^2*theta;
        end
        predict_cov = exp(-predict_cov/vari_num);

        % predict base fval

        Y_pred = fval_reg_pred_nomlz*beta+predict_cov'*gama;

        % predict variance
        u__ = fval_reg_nomlz'*inv_cov*predict_cov-fval_reg_pred_nomlz';
        Var_pred = sigma_sq*...
            (1+u__'/FTRF*u__+...
            -predict_cov'*inv_cov*predict_cov);

        % normalize data
        Y_pred = Y_pred*stdD_Y+aver_Y;
        Var_pred = diag(Var_pred)*stdD_Y*stdD_Y;
    end

    function F_reg = regZero(X)
        % zero order base funcion
        %
        F_reg = ones(size(X, 1), 1); % zero
    end

    function F_reg = regLinear(X)
        % first order base funcion
        %
        F_reg = [ones(size(X, 1), 1), X]; % linear
    end
end

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

[predict_function_fval,radialbasis_model_fval] = interpRadialBasisPreModel...
    (x_list,fval_list,basis_function);

if ~isempty(con_list)
    predict_function_con = cell(size(con_list,2),1);
    radialbasis_model_con = struct('X',[],'Y',[],...
        'radialbasis_matrix',[],'beta',[],...
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
        'radialbasis_matrix',[],[],'beta',[],...
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

[beta,rdibas_matrix] = interpRadialBasis...
    (X_dis,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function = @(X_predict) interpRadialBasisPredictor...
    (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,beta,basis_function);

radialbasis_model.X = X;
radialbasis_model.Y = Y;
radialbasis_model.radialbasis_matrix = rdibas_matrix;
radialbasis_model.beta = beta;

radialbasis_model.aver_X = aver_X;
radialbasis_model.stdD_X = stdD_X;
radialbasis_model.aver_Y = aver_Y;
radialbasis_model.stdD_Y = stdD_Y;
radialbasis_model.basis_function = basis_function;

radialbasis_model.predict_function = predict_function;

% abbreviation:
% num: number,pred: predict,vari: variable
    function [beta,rdibas_matrix] = interpRadialBasis...
            (X_dis,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix = rdibas_matrix+eye(x_number)*1e-6;

        % solve beta
        beta = rdibas_matrix\Y;
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
function [x_list,fval_list,con_list,coneq_list] = dataLibraryWrite...
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

function [x_list,fval_list,con_list,coneq_list] = dataLibraryRead...
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

function [data_library,x_list,fval_list,con_list,coneq_list,vio_list] = dataLibraryUpdata...
    (data_library,x_list,fval_list,con_list,coneq_list,vio_list)
% updata data to exist data library
%
x_list = [data_library.x_list;x_list];
data_library.x_list = x_list;
fval_list = [data_library.fval_list;fval_list];
data_library.fval_list = fval_list;
if ~isempty(data_library.con_list)
    con_list = [data_library.con_list;con_list];
    data_library.con_list = con_list;
end
if ~isempty(data_library.coneq_list)
    coneq_list = [data_library.coneq_list;coneq_list];
    data_library.coneq_list = coneq_list;
end
if ~isempty(data_library.vio_list)
    vio_list = [data_library.vio_list;vio_list];
    data_library.vio_list = vio_list;
end
end

function [x_list,fval_list,con_list,coneq_list,vio_list] = dataLibraryLoad...
    (data_library)
% updata data to exist data library
%
x_list = data_library.x_list;
fval_list = data_library.fval_list;
con_list = data_library.con_list;
coneq_list = data_library.coneq_list;
vio_list = data_library.vio_list;
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
