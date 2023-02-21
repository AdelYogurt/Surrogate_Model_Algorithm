clc;
clear;
close all hidden;

data_library_name='optimal_data_library';

benchmark=BenchmarkFunction();

variable_number=2;
object_function=@(x) benchmark.single2DObject(x);
object_function_LF=@(x) benchmark.single2DObjectLow(x);
A=[];
B=[];
Aeq=[];
Beq=[];
low_bou=[-5,-5];
up_bou=[5,5];
nonlcon_function=[];
nonlcon_function_LF=[];
cheapcon_function=[];

delete([data_library_name,'.txt']);
delete('result_total.txt');

[x_best,fval_best,NFE,output]=optimalSurrogateSADEKTS...
    (object_function,object_function_LF,variable_number,low_bou,up_bou,...
    nonlcon_function,nonlcon_function_LF,cheapcon_function,[],[],200,100);
result_x_best=output.result_x_best;
result_fval_best=output.result_fval_best;

figure(1);
plot(result_fval_best);

figure(2);
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);
scatter3(x_list(:,1),x_list(:,2),fval_list);

xlabel('X');
ylabel('Y');
zlabel('Z');

function [x_best,fval_best,NFE,output]=optimalSurrogateSADEKTS...
    (object_function,object_function_LF,variable_number,low_bou,up_bou,....
    nonlcon_function,nonlcon_function_LF,cheapcon_function,...
    model_function,model_function_LF,....
    NFE_max,iteration_max,torlance,nonlcon_torlance)
% SADE-KTS optimization algorithm
% adding knowledge-transfer-based sampling to original SADE algorithm 
%
% both nonlcon_function and cheapcon_function format is [con,coneq]
% model_function should output fval, format is [fval,con,coneq]
%
% referance: [1] LONG T, YE N, SHI R, et al. Surrogate-Assisted
% Differential Evolution Using Knowledge-Transfer-Based Sampling for
% Expensive Optimization Problems [J]. AIAA Journal, 2021, 60(1-16.
%
% Copyright 2022 Adel
%
if nargin < 14 || isempty(nonlcon_torlance)
    nonlcon_torlance=1e-3;
    if nargin < 13 || isempty(torlance)
        torlance=1e-3;
        if nargin < 12
            iteration_max=[];
            if nargin < 11
                NFE_max=[];
            end
        end
    end
end

if nargin < 10
    model_function_LF=[];
    if nargin < 9
        model_function=[];
        if nargin < 8
            cheapcon_function=[];
            if nargin < 7
                nonlcon_function_LF=[];
                if nargin < 6
                    nonlcon_function=[];
                end
            end
        end
    end
end

data_library_name='optimal_data_library';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

% hyper parameter
population_number=min(100,10*variable_number);
elite_rate=0.2;
correction_factor=0.5;

% generate knowledge transform data base
if isempty(model_function_LF)
    model_function_LF=@(x) modelFunction(x,object_function_LF,nonlcon_function_LF);
end
[~,x_list,~]=getLatinHypercube...
    (min(100,10*variable_number),variable_number,[],...
    low_bou,up_bou,cheapcon_function);
dataLibraryUpdata...
    (data_library_name,model_function_LF,x_list);

% Knowledge-Transfer-Based Sampling Method
x_initial_list=getInitialSample...
    (population_number,data_library_name,variable_number,...
    low_bou,up_bou,...
    elite_rate,correction_factor);

[x_best,fval_best,NFE,output]=optimalSurrogateSADE...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance,x_initial_list);

    function [fval,con,coneq]=modelFunction(x,object_function,nonlcon_function)
        % model function, concertrate fval, con, coneq into one function
        %
        if nargin < 3 || isempty(nonlcon_function)
            con=[];
            coneq=[];
        else
            [con,coneq]=nonlcon_function(x);
        end
        fval=object_function(x);
    end
end

function [x_best,fval_best,NFE,output]=optimalSurrogateSADE...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance,x_initial_list)
% KRG-CDE optimization algorithm
%
% referance: [1] Ҷ���, ����, �����, et al.
% ����Kriging����ģ�͵�Լ����ֽ����㷨 [J]. ����ѧ��, 2021, 42(6): 13.
%
% Copyright 2022 Adel
%
if nargin < 12
    x_initial_list=[];
    if nargin < 11 || isempty(nonlcon_torlance)
        nonlcon_torlance=1e-3;
        if nargin < 10 || isempty(torlance)
            torlance=1e-3;
            if nargin < 9
                iteration_max=[];
                if nargin < 8
                    NFE_max=[];
                end
            end
        end
    end
end

if nargin < 7
    model_function=[];
    if nargin < 6
        cheapcon_function=[];
        if nargin < 5
            nonlcon_function=[];
        end
    end
end

DRAW_FIGURE_FLAG=0; % whether draw data
INFORMATION_FLAG=1; % whether print data
CONVERGENCE_JUDGMENT_FLAG=0; % whether judgment convergence

% hyper parameter
population_number=min(100,10*variable_number);
RBF_number=max(100,(variable_number+1)*(variable_number+2)/2);

scaling_factor=0.8; % F
cross_rate=0.8;

protect_range=1e-4;

data_library_name='optimal_data_library';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;

% if do not input model_function, generate model_function
if isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end

iteration=iteration+1;

% step 2
% generate initial sample x_list
if isempty(x_initial_list)
    [~,x_updata_list,~]=getLatinHypercube...
        (population_number,variable_number,[],low_bou,up_bou,cheapcon_function);
else
    x_updata_list=x_initial_list;
end

% detech expensive constraints
if ~isempty(x_updata_list)
    [~,con,coneq]=dataLibraryUpdata...
        (data_library_name,model_function,x_updata_list(1,:));NFE=NFE+1;
    x_updata_list=x_updata_list(2:end,:);
else
    [~,con,coneq]=dataLibraryLoad(data_library_name,low_bou,up_bou);
end
if ~isempty(con) || ~isempty(coneq)
    expensive_nonlcon_flag=1;
else
    expensive_nonlcon_flag=0;
end

% NFE and iteration setting
if isempty(NFE_max)
    if expensive_nonlcon_flag
        NFE_max=50*variable_number;
    else
        NFE_max=20*variable_number;
    end
end

if isempty(iteration_max)
    if expensive_nonlcon_flag
        iteration_max=50*variable_number;
    else
        iteration_max=20*variable_number;
    end
end

result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);

% import data from data library
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);

% updata data library by x_list
[fval_updata_list,con_updata_list,coneq_updata_list]=dataLibraryUpdata...
    (data_library_name,model_function,x_updata_list);NFE=NFE+size(x_updata_list,1);
x_list=[x_list;x_updata_list];
fval_list=[fval_list;fval_updata_list];
con_list=[con_list;con_updata_list];
coneq_list=[coneq_list;coneq_updata_list];

search_flag=0; % global search or local search, 0 is global and 1 is local
while ~done
    infor_search_flag=search_flag;
    % nomalization con with average
    fval_max=mean(abs(fval_list),1);
    fval_nomlz_list=fval_list./fval_max*10;
    if ~isempty(con_list)
        con_max_list=mean(abs(con_list),1);
        con_nomlz_list=con_list./con_max_list*10;
    else
        con_max_list=[];
        con_nomlz_list=[];
    end
    if ~isempty(coneq_list)
        coneq_max_list=mean(abs(coneq_list),1);
        coneq_nomlz_list=coneq_list./coneq_max_list*10;
    else
        coneq_max_list=[];
        coneq_nomlz_list=[];
    end
    
    if search_flag == 0
        % global search
        [x_global_infill,feasiable_index_list]=searchGlobal...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
            population_number,scaling_factor,cross_rate,...
            expensive_nonlcon_flag,DRAW_FIGURE_FLAG);
        
        [x_global_infill,fval_global_infill,con_global_infill,coneq_global_infill,NFE_p]=dataLibraryUpdataProtect...
            (data_library_name,model_function,x_global_infill',...
            x_list,fval_list,con_list,coneq_list,...
            low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
        
        x_list=[x_list;x_global_infill];
        fval_list=[fval_list;fval_global_infill];
        con_list=[con_list;con_global_infill];
        coneq_list=[coneq_list;coneq_global_infill];
        
        % whether impove pupolation, if imporve, continue global
        search_flag=1;
        if isempty(feasiable_index_list)
            min_con=min(con_list,[],1);
            min_coneq=min(abs(coneq_list),[],1);
            if ~isempty(con_list)
                if ~sum(con_global_infill > min_con) % imporve
                    search_flag=0;
                end
            end
            if ~isempty(coneq_list)
                if ~sum(abs(coneq_global_infill) > min_coneq) % imporve
                    search_flag=0;
                end
            end
        else
            min_fval=min(fval_list,[],1);
            if fval_global_infill < min_fval  % imporve
                search_flag=0;
            end
        end
    else
        % local search
        x_local_infill=searchLocal...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
            population_number,RBF_number,...
            DRAW_FIGURE_FLAG);
        
        [x_local_infill,fval_local_infill,con_local_infill,coneq_local_infill,NFE_p]=dataLibraryUpdataProtect...
            (data_library_name,model_function,x_local_infill',...
            x_list,fval_list,con_list,coneq_list,...
            low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
        
        x_list=[x_list;x_local_infill];
        fval_list=[fval_list;fval_local_infill];
        con_list=[con_list;con_local_infill];
        coneq_list=[coneq_list;coneq_local_infill];
        
        % whether impove pupolation, if imporve, continue local
        search_flag=0;
        if isempty(feasiable_index_list)
            min_con=min(con_list,[],1);
            min_coneq=min(abs(coneq_list),[],1);
            if ~isempty(con_list)
                if ~sum(con_global_infill > min_con) % imporve
                    search_flag=1;
                end
            end
            if isempty(coneq_list)
                if ~sum(abs(coneq_global_infill) > min_coneq) % imporve
                    search_flag=1;
                end
            end
        else
            min_fval=min(fval_list,[],1);
            if fval_global_infill < min_fval  % imporve
                search_flag=1;
            end
        end
    end
    
    % find best result to record
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);
    
    if INFORMATION_FLAG
        fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
        if infor_search_flag == 0
            fprintf('global x:          %s\n',num2str(x_global_infill));
            fprintf('global value:      %f\n',fval_global_infill);
            fprintf('global violation:  %s  %s\n',num2str(con_global_infill),num2str(coneq_global_infill));
        else
            fprintf('local  x:          %s\n',num2str(x_local_infill));
            fprintf('local  value:      %f\n',fval_local_infill);
            fprintf('local  violation:  %s  %s\n',num2str(con_local_infill),num2str(coneq_local_infill));
        end
        fprintf('\n');
    end
    
    result_x_best(iteration,:)=x_best';
    result_fval_best(iteration,:)=fval_best;
    iteration=iteration+1;
    
    % forced interrupt
    if iteration > iteration_max || NFE >= NFE_max
        done=1;
    end
    
    % convergence judgment
    if CONVERGENCE_JUDGMENT_FLAG
        if (iteration > 1 && ...
                abs((fval_best-fval_best_old)/fval_best_old) < torlance)
            done=1;
            if ~isempty(con_best)
                if sum(con_best > nonlcon_torlance)
                    done=0;
                end
            end
            if ~isempty(coneq_best)
                if sum(abs(coneq_best) > nonlcon_torlance)
                    done=0;
                end
            end
        end
    end
    fval_best_old=fval_best;
    
end

result_x_best=result_x_best(1:iteration-1,:);
result_fval_best=result_fval_best(1:iteration-1);

output.result_x_best=result_x_best;
output.result_fval_best=result_fval_best;

    function [x_global_infill,feasiable_index_list]=searchGlobal...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
            population_number,scaling_factor,cross_rate,...
            expensive_nonlcon_flag,DRAW_FIGURE_FLAG)
        % find global infill point function
        %
        
        % step 4
        % updata kriging model and function
        [kriging_model_fval,kriging_model_con,kriging_model_coneq,output_kriging]=getKrigingModel...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list);
        object_function_surrogate=output_kriging.object_function_surrogate;
        object_function_variance=output_kriging.object_function_variance;
        nonlcon_function_surrogate=output_kriging.nonlcon_function_surrogate;
        nonlcon_function_variance=output_kriging.nonlcon_function_variance;
        
        % step 5
        % rank x_list data
        [x_rank_list,~,~,~]=rankData...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            cheapcon_function,nonlcon_torlance);
        
        % step 6
        % only the first population_number will be use
        x_best_popu_list=x_rank_list(1:population_number,:);
        
        % differ evolution mutations
        X_new_R1=differEvolutionRand...
            (low_bou,up_bou,x_best_popu_list,scaling_factor,population_number,1);
        X_new_R2=differEvolutionRand...
            (low_bou,up_bou,x_best_popu_list,scaling_factor,population_number,2);
        X_new_CR=differEvolutionCurrentRand...
            (low_bou,up_bou,x_best_popu_list,scaling_factor);
        X_new_CB=differEvolutionCurrentBest...
            (low_bou,up_bou,x_best_popu_list,scaling_factor,1);
        
        % differ evolution crossover
        X_new_R1=differEvolutionCrossover...
            (low_bou,up_bou,x_best_popu_list,X_new_R1,cross_rate);
        X_new_R2=differEvolutionCrossover...
            (low_bou,up_bou,x_best_popu_list,X_new_R2,cross_rate);
        X_new_CR=differEvolutionCrossover...
            (low_bou,up_bou,x_best_popu_list,X_new_CR,cross_rate);
        X_new_CB=differEvolutionCrossover...
            (low_bou,up_bou,x_best_popu_list,X_new_CB,cross_rate);
        
        % find global infill point base kriging model from offspring X
        x_DE_list=[X_new_R1;X_new_R2;X_new_CR;X_new_CB];
        fval_DE_list=zeros(4*population_number,1);
        fval_var_DE_list=zeros(4*population_number,1);
        if ~isempty(con_nomlz_list)
            con_DE_list=zeros(4*population_number,size(con_nomlz_list,2));
            con_var_DE_list=zeros(4*population_number,size(con_nomlz_list,2));
        end
        if ~isempty(coneq_nomlz_list)
            coneq_DE_list=zeros(4*population_number,size(coneq_nomlz_list,2));
            coneq_var_DE_list=zeros(4*population_number,size(coneq_nomlz_list,2));
        end
        
        % evaluate each x_offspring fval and constraints
        for x_index=1:4*population_number
            x=x_DE_list(x_index,:);
            fval_DE_list(x_index,1)=object_function_surrogate(x);
            fval_var_DE_list(x_index,1)=object_function_variance(x);
        end
        feasiable_index_list=[];
        if expensive_nonlcon_flag
            if ~isempty(nonlcon_function_surrogate)
                [con_temp,coneq_temp]=nonlcon_function_surrogate(x);
                if ~isempty(con_temp)
                    con_DE_list(x_index,:)=con_temp';
                end
                if ~isempty(coneq_temp)
                    coneq_DE_list(x_index,:)=coneq_temp';
                end
            end
            if ~isempty(nonlcon_function_variance)
                [con_var,coneq_var]=nonlcon_function_variance(x);
                if ~isempty(con_var)
                    con_var_DE_list(x_index,:)=con_var';
                end
                if ~isempty(coneq_var)
                    coneq_var_DE_list(x_index,:)=coneq_var';
                end
            end
            if max([con_temp;abs(coneq_temp)]) < nonlcon_torlance
                feasiable_index_list=[feasiable_index_list,x_index];
            end
        else
            feasiable_index_list=1:4*population_number;
        end
        
        % if have feasiable_index_list, only use feasiable to choose
        if ~isempty(feasiable_index_list)
            x_DE_list=x_DE_list(feasiable_index_list,:);
            fval_DE_list=fval_DE_list(feasiable_index_list);
            fval_var_DE_list=fval_var_DE_list(feasiable_index_list);
            if ~isempty(con_nomlz_list)
                con_DE_list=con_DE_list(feasiable_index_list,:);
                con_var_DE_list=con_var_DE_list(feasiable_index_list,:);
            end
            if ~isempty(coneq_nomlz_list)
                coneq_DE_list=coneq_DE_list(feasiable_index_list,:);
                coneq_var_DE_list=coneq_var_DE_list(feasiable_index_list,:);
            end
        end
        
        % base on constaints improve select global infill
        % lack process of equal constraints
        if isempty(feasiable_index_list)
            con_DE_base=max(min(con_DE_list,[],1),0);
            con_impove_probability_list=sum(...
                normpdf((con_DE_base-con_DE_list)./con_var_DE_list),2);
            [~,con_best_index]=max(con_impove_probability_list);
            con_best_index=con_best_index(1);
            x_global_infill=x_DE_list(con_best_index,:)';
        end
        
        % base on best fitness select global infill
        if ~isempty(feasiable_index_list)
            fval_DE_min=min(fval_DE_list,[],1);
            fval_DE_max=max(fval_DE_list,[],1);
            fval_var_DE_min=min(fval_var_DE_list,[],1);
            fval_var_DE_max=max(fval_var_DE_list,[],1);
            DE_fitness_list=-(fval_DE_list-fval_DE_min)/(fval_DE_max-fval_DE_min)+...
                0.5*(fval_var_DE_list-fval_var_DE_min)/(fval_var_DE_max-fval_var_DE_min);
            [~,fitness_best_index]=max(DE_fitness_list);
            fitness_best_index=fitness_best_index(1);
            x_global_infill=x_DE_list(fitness_best_index,:)';
        end
        
        if DRAW_FIGURE_FLAG && variable_number < 3
            line(x_global_infill(1),x_global_infill(2),'Marker','o','color','r')
            interpVisualize(kriging_model_fval,low_bou,up_bou);
        end
    end
    function x_local_infill=searchLocal...
            (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
            variable_number,low_bou,up_bou,cheapcon_function,nonlcon_torlance,...
            population_number,RBF_number,...
            DRAW_FIGURE_FLAG)
        % find local infill point function
        %
        
        % step 8
        % rand select initial local point from x_list
        x_index=randi(population_number);
        x_initial=x_list(x_index,:)';
        
        % rank x_list data
        %         [x_list_rank,~,~,~]=rankData...
        %             (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
        %             cheapcon_function,nonlcon_torlance);
        %         x_initial=x_list_rank(1,:)';
        
        % select nearest point to construct RBF
        RBF_number=min(RBF_number,size(x_list,1));
        distance=sum(((x_initial'-x_list)./(up_bou-low_bou)).^2,2);
        [~,index_list]=sort(distance);
        index_list=index_list(1:RBF_number);
        x_RBF_list=x_list(index_list,:);
        fval_RBF_nomlz_list=fval_nomlz_list(index_list,:);
        if ~isempty(con_nomlz_list)
            con_RBF_nomlz_list=con_nomlz_list(index_list,:);
        else
            con_RBF_nomlz_list=[];
        end
        if ~isempty(coneq_nomlz_list)
            coneq_RBF_nomlz_list=coneq_nomlz_list(index_list,:);
        else
            coneq_RBF_nomlz_list=[];
        end
        
        % get RBF model and function
        [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output_radialbasis]=getRadialBasisModel...
            (x_RBF_list,fval_RBF_nomlz_list,con_RBF_nomlz_list,coneq_RBF_nomlz_list);
        object_function_surrogate=output_radialbasis.object_function_surrogate;
        nonlcon_function_surrogate=output_radialbasis.nonlcon_function_surrogate;
        
        % get local infill point
        % obtian total constraint function
        if ~isempty(nonlcon_function_surrogate) || ~isempty(cheapcon_function)
            constraint_function=@(x) totalconFunction...
                (x,nonlcon_function_surrogate,cheapcon_function,nonlcon_torlance);
        else
            constraint_function=[];
        end
        fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp');
        x_local_infill=fmincon(object_function_surrogate,x_initial,[],[],[],[],...
            low_bou,up_bou,constraint_function,fmincon_options);
        
        if DRAW_FIGURE_FLAG && variable_number < 3
            line(x_local_infill(1),x_local_infill(2),'Marker','o','color','r')
            interpVisualize(radialbasis_model_fval,low_bou,up_bou);
        end
    end
    function [fval,con,coneq]=modelFunction(x,object_function,nonlcon_function)
        % model function, concertrate fval, con, coneq into one function
        %
        if nargin < 3 || isempty(nonlcon_function)
            con=[];
            coneq=[];
        else
            [con,coneq]=nonlcon_function(x);
        end
        fval=object_function(x);
    end
    function [con,coneq]=totalconFunction...
            (x,nonlcon_function,cheapcon_function,nonlcon_torlance)
        con=[];
        coneq=[];
        if ~isempty(nonlcon_function)
            [expencon,expenconeq]=nonlcon_function(x);
            con=[con;expencon-nonlcon_torlance];
            coneq=[coneq;expenconeq-nonlcon_torlance];
        end
        if ~isempty(cheapcon_function)
            [expencon,expenconeq]=cheapcon_function(x);
            con=[con;expencon];
            coneq=[coneq;expenconeq];
        end
    end
    function X_new=differEvolutionRand(low_bou,up_bou,X,F,x_number,rand_number)
        if nargin < 4
            rand_number=1;
            if nargin < 3
                x_number=1;
                if nargin < 2
                    error('differEvolutionRand: lack scaling factor F');
                end
            end
        end
        [x_number__,variable_number__]=size(X);
        X_new=zeros(x_number,variable_number__);
        for x_index__=1:x_number
            index__=randi(x_number__,2*rand_number+1,1);
            X_new(x_index__,:)=X(index__(1),:);
            for rand_index__=1:rand_number
                X_new(x_index__,:)=X_new(x_index__,:)+...
                    F*(X(index__(2*rand_index__),:)-X(index__(2*rand_index__+1),:));
                X_new(x_index__,:)=max(X_new(x_index__,:),low_bou);
                X_new(x_index__,:)=min(X_new(x_index__,:),up_bou);
            end
        end
    end
    function X_new=differEvolutionCurrentRand(low_bou,up_bou,X,F)
        [x_number__,variable_number__]=size(X);
        X_new=zeros(x_number__,variable_number__);
        for x_index__=1:x_number__
            index__=randi(x_number__,3,1);
            X_new(x_index__,:)=X(x_index__,:)+...
                F*(X(index__(1),:)-X(x_index__,:)+...
                X(index__(2),:)-X(index__(3),:));
            X_new(x_index__,:)=max(X_new(x_index__,:),low_bou);
            X_new(x_index__,:)=min(X_new(x_index__,:),up_bou);
        end
    end
    function X_new=differEvolutionCurrentBest(low_bou,up_bou,X,F,x_best_index)
        [x_number__,variable_number__]=size(X);
        X_new=zeros(x_number__,variable_number__);
        for x_index__=1:x_number__
            index__=randi(x_number__,2,1);
            X_new(x_index__,:)=X(x_index__,:)+...
                F*(X(x_best_index,:)-X(x_index__,:)+...
                X(index__(1),:)-X(index__(2),:));
            X_new(x_index__,:)=max(X_new(x_index__,:),low_bou);
            X_new(x_index__,:)=min(X_new(x_index__,:),up_bou);
        end
    end
    function X_new=differEvolutionCrossover(low_bou,up_bou,X,V,C_R)
        if size(X,1)~=size(V,1)
            error('differEvolutionOffspring: size incorrect');
        end
        [x_number__,variable_number__]=size(X);
        X_new=X;
        rand_number=rand(x_number__,variable_number__);
        index__=find(rand_number < C_R);
        X_new(index__)=V(index__);
        for x_index__=1:x_number__
            X_new(x_index__,:)=max(X_new(x_index__,:),low_bou);
            X_new(x_index__,:)=min(X_new(x_index__,:),up_bou);
        end
    end
    function [x_updata_list,fval_updata_list,con_updata_list,coneq_updata_list,NFE]=dataLibraryUpdataProtect...
            (data_library_name,model_function,x_add_list,...
            x_list,fval_list,con_list,coneq_list,...
            low_bou,up_bou,protect_range)
        % function updata data with same_point_avoid protect
        % return fval
        % all list is x_number x variable_number matrix
        %
        NFE=0;
        x_updata_list=[];fval_updata_list=[];con_updata_list=[];coneq_updata_list=[];
        for x_index__=1:size(x_add_list,1)
            x_updata__=x_add_list(x_index__,:);
            
            % check x_potential if exist in data library
            % if not, updata data libraray
            distance__=sum(((x_updata__-x_list)./(low_bou-up_bou)).^2,2);
            [~,min_index__]=min(distance__);
            if distance__(min_index__) < protect_range^2
                x_updata__=x_list(min_index__,:);
                fval_updata__=fval_list(min_index__,:);
                con_updata__=[];
                coneq_updata__=[];
                if ~isempty(con_list)
                    con_updata__=(con_list(min_index__,:));
                end
                if ~isempty(coneq_list)
                    coneq_updata__=(coneq_list(min_index__,:));
                end
            else
                [fval_updata__,con_updata__,coneq_updata__]=dataLibraryUpdata...
                    (data_library_name,model_function,x_updata__);NFE=NFE+1;
            end
            x_updata_list=[x_updata_list;x_updata__];
            fval_updata_list=[fval_updata_list;fval_updata__];
            con_updata_list=[con_updata_list;con_updata__];
            coneq_updata_list=[coneq_updata_list;coneq_updata__];
        end
    end
end


%% KTS
function [x_list]=getInitialSample...
    (population_number,data_library_name,variable_number,...
    low_bou,up_bou,...
    elite_rate,correction_factor)
% use SVM to correct latin hypercubic by exist data
% Knowledge-Transfer-Based Sampling Method
%

% generate initial latin hypercubic which will be corrected
[X_initial,~,~]=getLatinHypercube...
    (population_number,variable_number,[],...
    low_bou,up_bou);

% import data from data library to rank data
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);

% KTS
if ~isempty(x_list)
    % rank x_list data
    [X_input,~,~,~]=rankData...
        (x_list,fval_list,con_list,coneq_list);
    
    [x_number,variable_number]=size(x_list);
    
    N_elite=round(x_number*elite_rate);
    
    % generate SVM model
    Y=[ones(N_elite,1);-ones(x_number-N_elite,1)];
    [SVM_predict_function,~]=classifySupportVectorMachine...
        (X_input,Y,10,low_bou,up_bou);
    
    % get predict value by SVM
    Y=zeros(population_number,1);
    index_list=[];
    for x_index=1:population_number
        Y(x_index)=SVM_predict_function(X_initial(x_index,:)');
        if Y(x_index) > 0
            index_list=[index_list;x_index];
        end
    end
    
    while isempty(index_list)
        [X,X_new,distance_min_normalize]=getLatinHypercube...
            (population_number,variable_number,[],...
            low_bou,up_bou);
        % get predict value by SVM
        Y=zeros(population_number,1);
        index_list=[];
        for x_index=1:population_number
            Y(x_index)=SVM_predict_function(X(x_index,:)');
            if Y(x_index) > 0
                index_list=[index_list;x_index];
            end
        end
    end
    
    % move X to nearest X_superior
    X_superior=X_initial(index_list,:);
    X_inferior=X_initial;
    X_inferior(index_list,:)=[];
    for x_index=1:size(X_inferior,1)
        x=X_inferior(x_index,:);
        distance=sqrt(sum((x-X_superior).^2,2));
        [~,index]=min(distance);
        x_superior=X_superior(index,:); % nearest x_superior
        x=x+correction_factor*(x_superior-x);
        X_inferior(x_index,:)=x;
    end
    
    x_list=[X_inferior;X_superior];
else
    x_list=X_initial;
end
end

%% auxiliary function
function [x_best,fval_best,con_best,coneq_best]=findMinRaw...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance)
% find min fval in raw data
% x_list, rank is variable
% con_list, rank is con
% coneq_list, rank is coneq
% function will find min fval in con==0
% if there was not feasible x, find min consum
%
con_best=[];
coneq_best=[];
max_nonlcon_list=zeros(size(x_list,1),1);
max_cheapcon_list=zeros(size(x_list,1),1);
% process expendsive con
if ~isempty(con_list)
    max_nonlcon_list=max(con_list,[],2);
end
if ~isempty(coneq_list)
    max_nonlcon_list=max(abs(coneq_list),[],2);
end

% add cheap con
for x_index=1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x_list(x_index,:)');
        max_cheapcon_list(x_index)=max_cheapcon_list(x_index)+...
            sum(max(con,0))+sum(coneq.*coneq);
    end
end

con_judge_list=(max_nonlcon_list > nonlcon_torlance)+...
    (max_cheapcon_list > 0);
index=find(con_judge_list == 0);
if ~isempty(index)
    % feasible x
    x_list=x_list(index,:);
    fval_list=fval_list(index);
    if ~isempty(con_list)
        con_list=con_list(index,:);
    end
    if ~isempty(coneq_list)
        coneq_list=coneq_list(index,:);
    end
    
    % min fval
    [fval_best,index_best]=min(fval_list);
    x_best=x_list(index_best,:)';
    if ~isempty(con_list)
        con_best=con_list(index_best,:)';
    end
    if ~isempty(coneq_list)
        coneq_best=coneq_list(index_best,:)';
    end
else
    % min consum
    [~,index_best]=min(max_nonlcon_list);
    fval_best=fval_list(index_best);
    x_best=x_list(index_best,:)';
    if ~isempty(con_list)
        con_best=con_list(index_best,:)';
    end
    if ~isempty(coneq_list)
        coneq_best=coneq_list(index_best,:)';
    end
end
end

function [x_list,fval_list,con_list,coneq_list]=rankData...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance)
% rank data base on feasibility rule
% infeasible is rank by sum of constraint
%
if nargin < 6 || isempty(nonlcon_torlance)
    nonlcon_torlance=0;
end
if nargin < 5
    cheapcon_function=[];
end

[x_number,~]=size(x_list);
con_sum_list=zeros(x_number,1);
for x_index=1:x_number
    if ~isempty(con_list)
        con_sum_list(x_index)=con_sum_list(x_index)+sum(max(con_list(x_index,:)-nonlcon_torlance,0));
    end
    if ~isempty(coneq_list)
        con_sum_list(x_index)=con_sum_list(x_index)+sum(abs(coneq_list(x_index,:))-nonlcon_torlance);
    end
end

% add cheap con
for x_index=1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x_list(x_index,:)');
        con_sum_list(x_index)=con_sum_list(x_index)+...
            sum(max(con,0))+sum(max(abs(coneq),0));
    end
end

% get index_list to rank x by rule
index_list=1; % ranked index of x_list
if con_sum_list(1)==0
    infeasible_index=2;
else
    infeasible_index=1;
end
for x_index=2:x_number
    if con_sum_list(x_index)==0
        % local data
        local_index=1;
        while local_index <= (infeasible_index-1) &&...
                fval_list(x_index) > fval_list(index_list(local_index))
            local_index=local_index+1;
        end
        index_list=[index_list(1:local_index-1,:);...
            x_index;index_list(local_index:end,:)];
        infeasible_index=infeasible_index+1;
    else
        % local data
        local_index=infeasible_index;
        while local_index <= size(index_list,1) &&...
                con_sum_list(x_index) > con_sum_list(index_list(local_index))
            local_index=local_index+1;
        end
        index_list=[index_list(1:local_index-1,:);...
            x_index;index_list(local_index:end,:)];
    end
end

% rank by index_list
x_list=x_list(index_list,:);
fval_list=fval_list(index_list);
if ~isempty(con_list)
    con_list=con_list(index_list,:);
end
if ~isempty(coneq_list)
    coneq_list=coneq_list(index_list,:);
end
end

%% SVM
function [predict_function,SVM_model]=classifySupportVectorMachine...
    (X,class,C,low_bou,up_bou,kernel_function)
% generate support vector machine model version 0
% version 0 use fmincon to get alpha
% only support binary classification, 0 and 1
% X, Y is x_number x variable_number matrix
% C is penalty factor, default is empty
% kernel_function default is gauss kernal function
% kernel_function should be @(x1,x2) ...
%
if nargin < 6
    kernel_function=[];
    if nargin < 3
        C=[];
    end
end

[x_number,variable_number]=size(X);

% normalization data
if nargin < 5
    up_bou=max(X);
    if nargin < 4
        low_bou=min(X);
    end
end
X_nomlz=(X-low_bou)./(up_bou-low_bou);

% transfer class into y
Y=class;
for x_index=1:x_number
    if Y(x_index) == 0
        Y(x_index)=-1;
    end
end

% default kernal function
if isempty(kernel_function)
    sigma=-100*log(1/sqrt(x_number))/variable_number^2;
    kernel_function=@(x1,x2) exp(-((x1-x2)'*(x1-x2))*sigma);
end

% initialization kernal function process X_cov
X_cov=zeros(x_number);
for rank_index=1:x_number
    for colume_index=1:rank_index-1
        X_cov(rank_index,colume_index)=...
            X_cov(colume_index,rank_index);
    end
    for colume_index=rank_index:x_number
        X_cov(rank_index,colume_index)=...
            kernel_function(X_nomlz(rank_index,:)',X_nomlz(colume_index,:)');
    end
end

% min SVM object function to get alpha
object_function_SVM=@(alpha) -objectFunctionSVM(alpha,X_cov,Y);
alpha_initial=ones(x_number,1)*0.5;
low_bou_fmincon=0*ones(x_number,1);
if isempty(C) || C==0
    up_bou_fmincon=[];
else
    up_bou_fmincon=C*ones(x_number,1);
end
Aeq=Y';
fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp');
alpha=fmincon(object_function_SVM,alpha_initial,...
    [],[],Aeq,0,low_bou_fmincon,up_bou_fmincon,[],fmincon_options);

% obtain other paramter
w=sum(alpha.*Y.*X_nomlz);
index_list=find(alpha > 1e-6);

alpha_Y=alpha.*Y;
alpha_Y_cov=X_cov*alpha_Y;
b=sum(Y(index_list)-alpha_Y_cov(index_list))/length(index_list);

% generate predict function
predict_function=@(x) classifySupportVectorMachinePredictor...
    (x,X_nomlz,Y,alpha,b,low_bou,up_bou,kernel_function);

% output model
SVM_model.X=X;
SVM_model.class=class;
SVM_model.Y=Y;
SVM_model.X_nomlz=X_nomlz;
SVM_model.low_bou=low_bou;
SVM_model.up_bou=up_bou;
SVM_model.alpha=alpha;
SVM_model.w=w;
SVM_model.b=b;
SVM_model.kernel_function=kernel_function;
SVM_model.predict_function=predict_function;

    function fval=objectFunctionSVM(alpha,X_inner_product,Y)
        % support vector machine maximum object function
        %
        alpha=alpha(:);
        alpha_Y__=alpha.*Y;
        fval=sum(alpha)-alpha_Y__'*X_inner_product*alpha_Y__/2;
    end
    function [predict_class,predict_fval]=classifySupportVectorMachinePredictor...
            (x,X_nomlz,Y,alpha,b,low_bou,up_bou,kernel_function)
        % predict value of x is 1 or -1
        % x input is colume vector
        %
        x=x(:);
        
        x_number__=size(X_nomlz,1);
        x_nomlz=(x-low_bou')./(up_bou'-low_bou');
        X_inner_product__=zeros(x_number__,1);
        for x_index__=1:x_number__
            X_inner_product__(x_index__)=...
                kernel_function(X_nomlz(x_index__,:)',x_nomlz);
        end
        predict_fval=sum(alpha.*Y.*X_inner_product__)+b;
        predict_fval=1/(1+exp(-predict_fval));
        if predict_fval > 0.5
            predict_class=1;
        else
            predict_class=0;
        end
    end
end

%% surrogate model
function [kriging_model_fval,kriging_model_con,kriging_model_coneq,output]=getKrigingModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create kriging model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%

[predict_function_fval,kriging_model_fval]=interpKrigingPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    kriging_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'base_function_list',[],'fval_regression',[],'covariance',[],'inv_convariance',[],...
        'theta',[],'beta',[],'gama',[],'sigma_sq',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'predict_function',[]);
    kriging_model_con=repmat(kriging_model_con,[size(con_list,2),1]);
    predict_function_con=cell(size(con_list,2),1);
    for con_index=1:size(con_list,2)
        [predict_function_con{con_index},kriging_model_con(con_index)]=interpKrigingPreModel...
            (x_list,con_list(:,con_index));
    end
else
    kriging_model_con=[];
    predict_function_con=[];
end

if ~isempty(coneq_list)
    kriging_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'base_function_list',[],'fval_regression',[],'covariance',[],'inv_convariance',[],...
        'theta',[],'beta',[],'gama',[],'sigma_sq',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'predict_function',[]);
    kriging_model_coneq=repmat(kriging_model_coneq,[size(coneq_list,2),1]);
    predict_function_coneq=cell(size(con_list,2),1);
    for coneq_index=1:size(coneq_list,2)
        [predict_function_coneq{coneq_index},kriging_model_coneq(coneq_index)]=interpKrigingPreModel...
            (x_list,coneq_list(:,coneq_index),theta_coneq,base_function_list);
    end
else
    kriging_model_coneq=[];
    predict_function_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,predict_function_fval);
object_function_variance=@(predict_x) objectFunctionVariance(predict_x,predict_function_fval);
if isempty(predict_function_con) && isempty(predict_function_coneq)
    nonlcon_function_surrogate=[];
    nonlcon_function_variance=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,predict_function_con,predict_function_coneq);
    nonlcon_function_variance=@(predict_x) nonlconFunctionVariance(predict_x,predict_function_con,predict_function_coneq);
end

output.object_function_surrogate=object_function_surrogate;
output.object_function_variance=object_function_variance;
output.nonlcon_function_surrogate=nonlcon_function_surrogate;
output.nonlcon_function_variance=nonlcon_function_variance;
output.x_list=x_list;
output.fval_list=fval_list;
output.con_list=con_list;
output.coneq_list=coneq_list;


    function fval=objectFunctionSurrogate...
            (predict_x,predict_function)
        [fval,~]=predict_function(predict_x);
    end
    function fval_var=objectFunctionVariance...
            (predict_x,predict_function)
        [~,fval_var]=predict_function(predict_x);
    end
    function [con,coneq]=nonlconFunctionSurrogate...
            (predict_x,predict_function_con,predict_function_coneq)
        if isempty(predict_function_con)
            con=[];
        else
            con=zeros(length(predict_function_con),1);
            for con_index__=1:length(predict_function_con)
                [con(con_index__),~]=predict_function_con{con_index__}....
                    (predict_x);
            end
        end
        if isempty(predict_function_coneq)
            coneq=[];
        else
            coneq=zeros(length(predict_function_coneq),1);
            for coneq_index__=1:length(predict_function_con)
                [coneq(coneq_index__),~]=predict_function_coneq{coneq_index__}...
                    (predict_x);
            end
        end
    end
    function [con_var,coneq_var]=nonlconFunctionVariance...
            (predict_x,predict_function_con,predict_function_coneq)
        if isempty(predict_function_con)
            con_var=[];
        else
            con_var=zeros(length(predict_function_con),1);
            for con_index__=1:length(predict_function_con)
                [~,con_var(con_index__)]=predict_function_con{con_index__}....
                    (predict_x);
            end
        end
        if isempty(predict_function_coneq)
            coneq_var=[];
        else
            coneq_var=zeros(length(predict_function_coneq),1);
            for coneq_index__=1:length(predict_function_con)
                [~,coneq_var(coneq_index__)]=predict_function_coneq{coneq_index__}...
                    (predict_x);
            end
        end
    end
end

function [predict_function,kriging_model]=interpKrigingPreModel...
    (X,Y,hyp)
% version 7, nomalization method is grassian
% add multi x_predict input support
% prepare model, optimal theta and calculation parameter
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
% theta=exp(hyp)
%
% input initial data X, Y, which are real data
%
% output is a kriging model, include predict_function...
% X, Y, base_function_list
%
% Copyright 2023.2 Adel
%
[x_number,variable_number]=size(X);
if nargin < 3
    hyp=zeros(1,variable_number);
end

% normalize data
aver_X=mean(X);
stdD_X=std(X);
aver_Y=mean(Y);
stdD_Y=std(Y);
index__=find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__)=1; end
index__=find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__)=1; end
X_nomlz=(X-aver_X)./stdD_X;
Y_nomlz=(Y-aver_Y)./stdD_Y;

% initial X_dis_sq
X_dis_sq=zeros(x_number,x_number,variable_number);
for variable_index=1:variable_number
    X_dis_sq(:,:,variable_index)=...
        (X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end

% regression function define
% notice reg_function process no normalization data
reg_function=@(X) regZero(X);
% reg_function=@(X) regLinear(X);

% calculate reg
fval_reg_nomlz=(reg_function(X)-0)./1;

% optimal to get hyperparameter
fmincon_option=optimoptions('fmincon','Display','none',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',false);
low_bou_hyp=-4*ones(1,variable_number);
up_bou_hyp=4*ones(1,variable_number);
object_function_hyp=@(hyp) objectNLLKriging...
    (X_dis_sq,Y_nomlz,x_number,variable_number,hyp,fval_reg_nomlz);

% [fval,gradient]=object_function_hyp(hyp)
% [~,gradient_differ]=differ(object_function_hyp,hyp)

% drawFunction(object_function_hyp,low_bou_hyp,up_bou_hyp);

hyp=fmincon...
    (object_function_hyp,hyp,[],[],[],[],low_bou_hyp,up_bou_hyp,[],fmincon_option);

% get parameter
[covariance,inv_covariance,beta,sigma_sq]=interpKriging...
    (X_dis_sq,Y_nomlz,x_number,variable_number,exp(hyp),fval_reg_nomlz);
gama=inv_covariance*(Y_nomlz-fval_reg_nomlz*beta);
FTRF=fval_reg_nomlz'*inv_covariance*fval_reg_nomlz;

% initialization predict function
predict_function=@(X_predict) interpKrigingPredictor...
    (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,exp(hyp),beta,gama,sigma_sq,...
    inv_covariance,fval_reg_nomlz,FTRF,reg_function);

kriging_model.X=X;
kriging_model.Y=Y;
kriging_model.X_normalize=X_nomlz;
kriging_model.Y_normalize=Y_nomlz;
kriging_model.fval_regression=fval_reg_nomlz;
kriging_model.covariance=covariance;
kriging_model.inv_covariance=inv_covariance;

kriging_model.hyp=hyp;
kriging_model.beta=beta;
kriging_model.gama=gama;
kriging_model.sigma_sq=sigma_sq;
kriging_model.aver_X=aver_X;
kriging_model.stdD_X=stdD_X;
kriging_model.aver_Y=aver_Y;
kriging_model.stdD_Y=stdD_Y;

kriging_model.predict_function=predict_function;

% abbreviation:
% num: number, pred: predict, vari: variable, hyp: hyper parameter
% NLL: negative log likelihood
    function [fval,gradient]=objectNLLKriging...
            (X_dis_sq,Y,x_num,vari_num,hyp,F_reg)
        % function to minimize sigma_sq
        %
        theta=exp(hyp);
        [cov,inv_cov,~,sigma2,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg);

        % calculation negative log likelihood
        L=chol(cov)';
        fval=x_num/2*log(sigma2)+sum(log(diag(L)));

        % calculate gradient
        if nargout > 1
            % gradient
            gradient=zeros(vari_num,1);
            for vari_index=1:vari_num
                dcov_dtheta=-(X_dis_sq(:,:,vari_index).*cov)*theta(vari_index)/vari_num;

                dinv_cov_dtheta=...
                    -inv_cov*dcov_dtheta*inv_cov;

                dinv_FTRF_dtheta=-inv_FTRF*...
                    (F_reg'*dinv_cov_dtheta*F_reg)*...
                    inv_FTRF;
                
                dmiu_dtheta=dinv_FTRF_dtheta*(F_reg'*inv_cov*Y)+...
                    inv_FTRF*(F_reg'*dinv_cov_dtheta*Y);
                
                dY_Fmiu_dtheta=-F_reg*dmiu_dtheta;

                dsigma2_dtheta=(dY_Fmiu_dtheta'*inv_cov*Y_Fmiu+...
                    Y_Fmiu'*dinv_cov_dtheta*Y_Fmiu+...
                    Y_Fmiu'*inv_cov*dY_Fmiu_dtheta)/x_num;
                
                dlnsigma2_dtheta=1/sigma2*dsigma2_dtheta;

                dlndetR=trace(inv_cov*dcov_dtheta);

                gradient(vari_index)=x_num/2*dlnsigma2_dtheta+0.5*dlndetR;
            end
        end
    end
    function [cov,inv_cov,beta,sigma_sq,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg)
        % kriging interpolation kernel function
        % Y(x)=beta+Z(x)
        %
        cov=zeros(x_num,x_num);
        for vari_index=1:vari_num
            cov=cov+X_dis_sq(:,:,vari_index)*theta(vari_index);
        end
        cov=exp(-cov/vari_num)+eye(x_num)*1e-3;

        % coefficient calculation
        inv_cov=cov\eye(x_num);
        inv_FTRF=(F_reg'*inv_cov*F_reg)\eye(size(F_reg,2));

        % basical bias
        beta=inv_FTRF*(F_reg'*inv_cov*Y);
        Y_Fmiu=Y-F_reg*beta;
        sigma_sq=(Y_Fmiu'*inv_cov*Y_Fmiu)/x_num;
        
    end
    function [Y_pred,Var_pred]=interpKrigingPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,theta,beta,gama,sigma_sq,...
            inv_cov,fval_reg_nomlz,FTRF,reg_function)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num,~]=size(X_pred);
        fval_reg_pred=reg_function(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;
        fval_reg_pred_nomlz=(fval_reg_pred-0)./1;
        
        % predict covariance
        predict_cov=zeros(x_num,x_pred_num);
        for vari_index=1:vari_num
            predict_cov=predict_cov+...
                (X_nomlz(:,vari_index)-X_pred_nomlz(:,vari_index)').^2*theta(vari_index);
        end
        predict_cov=exp(-predict_cov/vari_num);

        % predict base fval
        
        Y_pred=fval_reg_pred_nomlz*beta+predict_cov'*gama;
        
        % predict variance
        u__=fval_reg_nomlz'*inv_cov*predict_cov-fval_reg_pred_nomlz';
        Var_pred=sigma_sq*...
            (1+u__'/FTRF*u__+...
            -predict_cov'*inv_cov*predict_cov);
        
        % normalize data
        Y_pred=Y_pred*stdD_Y+aver_Y;
        Var_pred=diag(Var_pred)*stdD_Y*stdD_Y;
    end
    function F_reg=regZero(X)
        % zero order base funcion
        %
        F_reg=ones(size(X,1),1); % zero
    end
    function F_reg=regLinear(X)
        % first order base funcion
        %
        F_reg=[ones(size(X,1),1),X]; % linear
    end
end

function [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output]=getRadialBasisModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create radialbasis model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
radialbasis_model_fval=interpRadialBasisPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    radialbasis_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'radialbasis_matrix',[],'inv_radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_con=repmat(radialbasis_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        radialbasis_model_con(con_index)=interpRadialBasisPreModel...
            (x_list,con_list(:,con_index));
    end
else
    radialbasis_model_con=[];
end

if ~isempty(coneq_list)
    radialbasis_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'radialbasis_matrix',[],'inv_radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_coneq=repmat(radialbasis_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        radialbasis_model_coneq(coneq_index)=interpRadialBasisPreModel...
            (x_list,coneq_list(:,coneq_index));
    end
else
    radialbasis_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,radialbasis_model_fval);
if isempty(radialbasis_model_con) && isempty(radialbasis_model_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,radialbasis_model_con,radialbasis_model_coneq);
end

output.object_function_surrogate=object_function_surrogate;
output.nonlcon_function_surrogate=nonlcon_function_surrogate;
output.x_list=x_list;
output.fval_list=fval_list;
output.con_list=con_list;
output.coneq_list=coneq_list;

    function fval=objectFunctionSurrogate...
            (predict_x,ERBF_model_fval)
        fval=ERBF_model_fval.predict_function(predict_x);
    end
    function [con,coneq]=nonlconFunctionSurrogate...
            (predict_x,ERBF_model_con,ERBF_model_coneq)
        if isempty(ERBF_model_con)
            con=[];
        else
            con=zeros(length(ERBF_model_con),1);
            for con_index__=1:length(ERBF_model_con)
                con(con_index__)=ERBF_model_con...
                    (con_index__).predict_function(predict_x);
            end
        end
        if isempty(ERBF_model_coneq)
            coneq=[];
        else
            coneq=zeros(length(ERBF_model_coneq),1);
            for coneq_index__=1:length(ERBF_model_con)
                coneq(coneq_index__)=ERBF_model_coneq...
                    (coneq_index__).predict_function(predict_x);
            end
        end
    end
end

function radialbasis_model=interpRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interp pre model function version 1
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
% beta is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
if nargin < 9
    X_nomlz=[];
    if nargin < 3
        basis_function=[];
    end
end

[x_number,variable_number]=size(X);

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__ = find(stdD_X == 0);
if ~isempty(index__),stdD_X(index__)=1;end
index__ = find(stdD_Y == 0);
if ~isempty(index__),stdD_Y(index__)=1;end
X_nomlz = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_nomlz = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

if isempty(basis_function)
    c=(prod(max(X_nomlz)-min(Y_nomlz))/x_number)^(1/variable_number);
    basis_function=@(r) exp(-(r'*r)/c);
    %     basis_function=@(r) sqrt(r'*r+c*c);
end

[beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
    (X_nomlz,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function=@(predict_x) interpRadialBasisPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    beta,basis_function,predict_x);

radialbasis_model.X=X;
radialbasis_model.Y=Y;
radialbasis_model.X_normalize=X_nomlz;
radialbasis_model.Y_normalize=Y_nomlz;
radialbasis_model.radialbasis_matrix=rdibas_matrix;
radialbasis_model.inv_radialbasis_matrix=inv_rdibas_matrix;

radialbasis_model.beta=beta;
radialbasis_model.aver_X=aver_X;
radialbasis_model.stdD_X=stdD_X;
radialbasis_model.aver_Y=aver_Y;
radialbasis_model.stdD_Y=stdD_Y;
radialbasis_model.basis_function=basis_function;

radialbasis_model.predict_function=predict_function;

    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
            (X,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix=zeros(x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                rdibas_matrix(rank_index__,colume_index__)=...
                    rdibas_matrix(colume_index__,rank_index__);
            end
            for colume_index__=rank_index__:x_number
                rdibas_matrix(rank_index__,colume_index__)=...
                    basis_function(X(rank_index__,:)'-X(colume_index__,:)');
            end
        end
        
        % stabilize matrix
        rdibas_matrix=rdibas_matrix+eye(x_number)*1e-6;
        
        inv_rdibas_matrix=inv(rdibas_matrix);
        beta=inv_rdibas_matrix*Y;
    end
    function [predict_y]=interpRadialBasisPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            beta,basis_function,predict_x)
        % radial basis function interp predict function
        % input predict_x and radialbasis_model model
        % predict_x is row vector
        % output the predict value
        %
        % Copyright 2022 Adel
        %
        if size(predict_x,2) > 1
            predict_x=predict_x';
        end
        
        [x_number__,~]=size(X_nomlz);
        
        % normalize data
        predict_x=(predict_x-aver_X')./stdD_X';
        
        % predict value
        X_inner_product__=zeros(x_number__,1);
        for index_i=1:x_number__
            X_inner_product__(index_i,1)=...
                basis_function(X_nomlz(index_i,:)'-predict_x);
        end
        
        % predict variance
        predict_y=beta'*X_inner_product__;
        
        % normalize data
        predict_y=predict_y*stdD_Y+aver_Y;
    end
end

%% data library
function [fval_list,con_list,coneq_list]=dataLibraryUpdata...
    (data_library_name,model_function,x_list)
% updata data library
%
fval_list=[];
con_list=[];
coneq_list=[];
[x_number,variable_number]=size(x_list);

if ~strcmp(data_library_name(end-3:end),'.txt')
    data_library_name=[data_library_name,'.txt'];
end

% store format
x_format_base='%.8e ';
fval_format_base='%.8e ';
x_format=repmat(x_format_base,1,variable_number);
file_optimalSurrogate_output = fopen(data_library_name,'a');
file_result = fopen('result_total.txt','a');

% updata
for x_index=1:x_number
    x=x_list(x_index,:)';
    [fval,con,coneq]=model_function(x);
    fval_list=[fval_list;fval];
    con_list=[con_list;con'];
    coneq_list=[coneq_list;coneq'];
    
    % write data to txt_optimalSurrogateSADEKTS
    fprintf(file_optimalSurrogate_output,'%d ',variable_number);
    fprintf(file_optimalSurrogate_output,'%d ',length(con));
    fprintf(file_optimalSurrogate_output,'%d ',length(coneq));
    fprintf(file_optimalSurrogate_output,x_format,x);
    fprintf(file_optimalSurrogate_output,fval_format_base,fval);
    fval_format=repmat(fval_format_base,1,length(con));
    fprintf(file_optimalSurrogate_output,fval_format,con);
    fval_format=repmat(fval_format_base,1,length(coneq));
    fprintf(file_optimalSurrogate_output,fval_format,coneq);
    fprintf(file_optimalSurrogate_output,'\n');
    
    % write data to txt_result
    fprintf(file_result,'%d ',variable_number);
    fprintf(file_result,'%d ',length(con));
    fprintf(file_result,'%d ',length(coneq));
    fprintf(file_result,x_format,x);
    fprintf(file_result,fval_format_base,fval);
    fval_format=repmat(fval_format_base,1,length(con));
    fprintf(file_result,fval_format,con);
    fval_format=repmat(fval_format_base,1,length(coneq));
    fprintf(file_result,fval_format,coneq);
    fprintf(file_result,'\n');
end

fclose(file_optimalSurrogate_output);
clear('file_optimalSurrogate_output');
fclose(file_result);
clear('file_result');
end

function [x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou)
% load data from data library
% low_bou, up_bou is range of data
%
if nargin < 3
    up_bou=inf;
    if nargin < 2
        low_bou=-inf;
        if nargin < 1
            error('dataLibraryLoad: lack data_library_name');
        end
    end
end

if ~strcmp(data_library_name(end-3:end),'.txt')
    data_library_name=[data_library_name,'.txt'];
end

if exist(data_library_name,'file')==2
    data_list=importdata(data_library_name);
    if ~isempty(data_list)
        % search whether exist point
        x_list=[];
        fval_list=[];
        con_list=[];
        coneq_list=[];
        
        for data_index=1:size(data_list,1)
            data=data_list(data_index,:);
            variable_number=data(1);
            con_number=data(2);
            coneq_number=data(3);
            
            x=data(4:3+variable_number);
            judge=sum(x < low_bou)+sum(x > up_bou);
            if ~judge
                x_list=[x_list;x];
                fval_list=[fval_list;data(4+variable_number)];
                con=data(5+variable_number:4+variable_number+con_number);
                if ~isempty(con)
                    con_list=[con_list;con];
                end
                coneq=data(6+variable_number+con_number:5+variable_number+con_number+coneq_number);
                if ~isempty(coneq)
                    coneq_list=[coneq_list;coneq];
                end
            end
        end
    else
        x_list=[];
        fval_list=[];
        con_list=[];
        coneq_list=[];
    end
else
    x_list=[];
    fval_list=[];
    con_list=[];
    coneq_list=[];
end
end

%% LHD
function [X,X_new,distance_min_nomlz]=getLatinHypercube...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou,cheapcon_function)
% generate sample sequence latin hypercube
% more uniform point distribution by simulating particle motion
% sample number is total point in area
% default low_bou is 0, up_bou is 1, cheapcon_function is []
% low_bou and up_bou is colume vector
% x in x_exist_list, x_list, supply_x_list is row vector
% x_exist_list should meet bou
%
% Copyright 2022 Adel
%
if nargin < 6
    cheapcon_function=[];
    if nargin < 5
        if nargin < 3
            X_exist=[];
            if nargin < 2
                error('getLatinHypercube: lack variable_number');
            end
        end
        low_bou=zeros(1,variable_number);
        up_bou=ones(1,variable_number);
    end
end
iteration_max=100;

% check x_exist_list if meet boundary
if ~isempty(X_exist)
    index=find(X_exist < low_bou);
    index=[index,find(X_exist > up_bou)];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list range error');
    end
    if size(X_exist,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
    X_exist_nomlz=(X_exist-low_bou)./(up_bou-low_bou);
else
    X_exist_nomlz=[];
end

% check input
if sample_number < 0
    X=[];
    X_new=[];
    distance_min_nomlz=[];
    return;
end

% check x_new_number
x_new_number=sample_number-size(X_exist,1);
if x_new_number < 0
    X=X_exist;
    X_new=[];
    distance_min_nomlz=getMinDistance(X_exist_nomlz);
    return;
end

low_bou_nomlz=zeros(1,variable_number);
up_bou_nomlz=ones(1,variable_number);

% get initial X_new_nomalize by lhsdesign
X_new_nomlz=rand(x_new_number,variable_number);
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);

% x is nomalize, so constraint function should change
if ~isempty(cheapcon_function)
    cheapcon_function=@(x) ...
        max(max(sample_number*cheapcon_function(x.*(up_bou-low_bou)+low_bou)+1,0),[],1);
end

iteration=0;
fval_list=zeros(x_new_number,1);
gradient_list=zeros(x_new_number,variable_number);
while iteration < iteration_max
    % change each x place by newton methods
    for x_index=1:x_new_number    
        % get gradient
        [fval_list(x_index,1),gradient_list(x_index,:)]=objectFunctionXPlace...
            (X_new_nomlz(x_index,:),[X_new_nomlz(1:x_index-1,:);X_new_nomlz(x_index+1:end,:);X_exist_nomlz],...
            sample_number,variable_number,low_bou_nomlz-0.1/variable_number,up_bou_nomlz+0.1/variable_number,cheapcon_function);
    end
    
    % normalize fval
    fval_list=fval_list/max(fval_list);
    for x_index=1:x_new_number
        C=fval_list(x_index,1)*distance_min_nomlz*(1-iteration/iteration_max);
        x=X_new_nomlz(x_index,:)+...
            -gradient_list(x_index,:)/...
            norm(gradient_list(x_index,:))*C;
        x=min(x,up_bou_nomlz);
        x=max(x,low_bou_nomlz);
        X_new_nomlz(x_index,:)=x;
    end
    
    iteration=iteration+1;
end
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);
X_new=X_new_nomlz.*(up_bou-low_bou)+low_bou;
X=[X_new;X_exist];

    function [fval,gradient]=objectFunctionXPlace...
            (x,X_surplus,sample_number,variable_number,low_bou,up_bou,cheapcon_function)
        % function describe distance between X and X_supply
        % x is colume vector and X_surplus is matrix which is num-1 x var
        % low_bou_limit__ and up_bou_limit__ is colume vector
        % variable in colume
        %
        a__=10/variable_number;
        a_bou__=30/sample_number;
        
        sign__=((x > X_surplus)-0.5)*2;
        
        xi__=-a__*(x-X_surplus).*sign__;
        sum_xi__=sum(xi__,2);
        psi__=a__*(low_bou-x)*a_bou__;
        zeta__=a__*(x-up_bou)*a_bou__;
        
%         exp_xi__=exp(xi__);
        exp_sum_xi__=exp(sum_xi__);
        exp_psi__=exp(psi__);
        exp_zeta__=exp(zeta__);
        
        % get fval
        fval=sum(exp_sum_xi__,1)+...
            sum(exp_psi__+exp_zeta__,2);
        
        % get gradient
        gradient=sum(-a__*sign__.*exp_sum_xi__,1)+...
            -a__*exp_psi__*a_bou__+...
            a__*exp_zeta__*a_bou__;
        
        if ~isempty(cheapcon_function)
            fval_con=cheapcon_function(x);
            fval=fval+fval_con;
            [gradient_con]=differ...
                (cheapcon_function,x,fval_con,variable_number);
            gradient=gradient+gradient_con;
        end
        
        function [gradient]=differ(differ_function,x,fval,variable_number,step)
            % differ function to get gradient and hessian
            % gradient is rank vector
            %
            if nargin < 5
                step=1e-6;
            end
            fval__=zeros(variable_number,2); % backward is 1, forward is 2
            gradient=zeros(1,variable_number);
            
            % fval and gradient
            for variable_index__=1:variable_number
                x_forward__=x;
                x_backward__=x;
                x_backward__(variable_index__)=x_backward__(variable_index__)-step;
                fval__(variable_index__,1)=differ_function(x_backward__);
                
                x_forward__(variable_index__)=x_forward__(variable_index__)+step;
                fval__(variable_index__,2)=differ_function(x_forward__);
                
                gradient(variable_index__)=...
                    (fval__(variable_index__,2)-fval__(variable_index__,1))/2/step;
            end
        end
    end
    function distance_min__=getMinDistance(x_list__)
        % get distance min from x_list
        %
        
        % sort x_supply_list_initial to decrese distance calculate times
        x_list__=sortrows(x_list__,1);
        sample_number__=size(x_list__,1);
        variable_number__=size(x_list__,2);
        distance_min__=variable_number__;
        for x_index__=1:sample_number__
            x_curr__=x_list__(x_index__,:);
            x_next_index__=x_index__ + 1;
            % first dimension only search in min_distance
            search_range__=variable_number__;
            while x_next_index__ <= sample_number__ &&...
                    (x_list__(x_next_index__,1)-x_list__(x_index__,1))^2 ...
                    < search_range__
                x_next__=x_list__(x_next_index__,:);
                distance_temp__=sum((x_next__-x_curr__).^2);
                if distance_temp__ < distance_min__
                    distance_min__ = distance_temp__;
                end
                if distance_temp__ < search_range__
                    search_range__ = distance_temp__;
                end
                x_next_index__=x_next_index__+1;
            end
        end
        distance_min__=sqrt(distance_min__);
    end
end