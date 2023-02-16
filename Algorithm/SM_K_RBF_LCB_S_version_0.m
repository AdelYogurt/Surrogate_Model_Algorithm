clc;
clear;
close all hidden;

data_library_name='optimalSurrogate_K_LCB_RBF_S_result';

% variable_number=2;
% object_function=@(x) functionGPObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=-3*ones(variable_number,1);
% up_bou=3*ones(variable_number,1);
% nonlcon_function=[];

variable_number=2;
object_function=@(x) functionPKObject(x);
A=[];
B=[];
Aeq=[];
Beq=[];
low_bou=-3*ones(variable_number,1);
up_bou=3*ones(variable_number,1);
nonlcon_function=[];
cheapcon_function=[];
model_function=[];

% variable_number=6;
% object_function=@(x) functionHNObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=zeros(variable_number,1);
% up_bou=ones(variable_number,1);
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

delete([data_library_name,'.txt']);
delete('result_total.txt');
[x_best,fval_best,NFE,output]=optimalSurrogateKLCBRBFS...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function)
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

function [x_best,fval_best,NFE,output]=optimalSurrogateKLCBRBFS...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance)
% surrogate base optimal KRG+LCB+RBF version 0
% use PCB function to decide which point was select to updata
% use kriging model to decrease model_function call
% and ues ga to optimal
% iteration the model an optimal
% x_best is colume vector by ga
% all function exchange x should be colume vector
% model_function should output fval, format is [fval,con,coneq]
% con or coneq can be colume vector if there was more than one constrain
%
% Copyright 2022 Adel
%
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

if nargin < 8
    NFE_max=[];
    if nargin < 7
        cheapcon_function=[];
        if nargin < 6
            nonlcon_function=[];
        end
    end
end
DRAW_FIGURE_FLAG=0; % whether draw data
INFORMATION_FLAG=1; % whether print data
CONVERGENCE_JUDGMENT_FLAG=0; % whether judgment convergence

sample_number_initial=5*variable_number;
RBF_number=6;

w_list=[1];

protect_range=1e-2;

data_library_name='optimalSurrogate_K_LCB_RBF_S_result';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;

% if do not input model_function, generate model_function
if nargin < 11 || isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end

iteration=iteration+1;

% check whether exist data
[x_list,~,~,~]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);

% generate latin hypercube sequence
[~,x_updata_list,~]=getLatinHypercube...
    (sample_number_initial,variable_number,x_list,...
    low_bou,up_bou,cheapcon_function);

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

% NFE_max and iteration setting
if expensive_nonlcon_flag
    if isempty(NFE_max)
        NFE_max=25*variable_number^2;
    end
    if isempty(iteration_max)
        iteration_max=100;
    end
    range_min=0.001;
else
    if isempty(NFE_max)
        NFE_max=10*variable_number^2;
    end
    if isempty(iteration_max)
        iteration_max=50;
    end
    range_min=0.01;
end

result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);

% import data from data library
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);

% updata X_updata into data library
[fval_updata_list,con_updata_list,coneq_updata_list]=dataLibraryUpdata...
    (data_library_name,model_function,x_updata_list);NFE=NFE+size(x_updata_list,1);

x_list=[x_list;x_updata_list];
fval_list=[fval_list;fval_updata_list];
con_list=[con_list;con_updata_list];
coneq_list=[coneq_list;coneq_updata_list];

while ~done
    % nomalization con average
    fval_max=mean(abs(fval_list),1);
    fval_nomlz_list=fval_list./fval_max;
    if ~isempty(con_list)
        con_max_list=mean(abs(con_list),1);
        con_nomlz_list=con_list./con_max_list;
    else
        con_max_list=[];
        con_nomlz_list=[];
    end
    if ~isempty(coneq_list)
        coneq_max_list=mean(abs(coneq_list),1);
        coneq_nomlz_list=coneq_list./coneq_max_list;
    else
        coneq_max_list=[];
        coneq_nomlz_list=[];
    end
    
    % get kriging model and function
    [kriging_model_fval,kriging_model_con,kriging_model_coneq,output_kriging]=getKrigingModel...
        (x_list,fval_list,con_list,coneq_list);
    object_function_surrogate=output_kriging.object_function_surrogate;
    object_function_variance=output_kriging.object_function_variance;
    nonlcon_function_surrogate=output_kriging.nonlcon_function_surrogate;
    nonlcon_function_variance=output_kriging.nonlcon_function_variance;
    
    % LCB guideline to obtain x_global_infill
    w=w_list(mod(iteration-1,length(w_list))+1);
    [x_global_infill,~,lcb_function]=findMinLCB...
        (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function_surrogate,cheapcon_function,...
        object_function_variance,w);
    
    % add x_global_infill
    [x_global_infill,fval_global_infill,con_global_infill,coneq_global_infill,NFE_p]=dataLibraryUpdataProtect...
        (data_library_name,model_function,x_global_infill',...
        x_list,fval_list,con_list,coneq_list,...
        low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
    
    x_list=[x_list;x_global_infill];
    fval_list=[fval_list;fval_global_infill];
    if ~isempty(con_list)
        con_list=[con_list;con_global_infill];
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_global_infill];
    end
    
    fval_nomlz_list=[fval_nomlz_list;fval_global_infill/fval_max];
    if ~isempty(con_nomlz_list)
        con_nomlz_list=[con_nomlz_list;con_global_infill./con_max_list];
    end
    if ~isempty(con_nomlz_list)
        coneq_nomlz_list=[coneq_nomlz_list;coneq_global_infill./coneq_max_list];
    end
    
    if DRAW_FIGURE_FLAG && variable_number < 3
        interpolationVisualize(kriging_model_fval,low_bou,up_bou);
        line(x_global_infill(1),x_global_infill(2),fval_global_infill,'Marker','o','color','r')
    end
    
    % rank x_list data
    [x_list_rank,~,~,~]=rankData...
        (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
        cheapcon_function,nonlcon_torlance);
    x_RBF_center=x_list_rank(1,:)';
    
    distance=sum(((x_RBF_center'-x_list)./(up_bou'-low_bou')).^2,2);
    [~,index_list]=sort(distance);
    index_list=index_list(1:RBF_number);
    x_RBF_list=x_list(index_list,:);
    fval_RBF_list=fval_list(index_list,:);
    if ~isempty(con_nomlz_list)
        con_RBF_list=con_list(index_list,:);
    else
        con_RBF_list=[];
    end
    if ~isempty(coneq_nomlz_list)
        coneq_RBF_list=coneq_list(index_list,:);
    else
        coneq_RBF_list=[];
    end
    
    % get RBF model and function
    %     [RadialBasis_model_fval,RadialBasis_model_con,RadialBasis_model_coneq,output_RadialBasis]=getRadialBasisModel...
    %         (x_RBF_list,fval_RBF_list,con_RBF_list,coneq_RBF_list);
    %     object_function_surrogate=output_RadialBasis.object_function_surrogate;
    %     nonlcon_function_surrogate=output_RadialBasis.nonlcon_function_surrogate;
    
    % get RPS model and function
    [respsurf_model_fval,respsurf_model_con,respsurf_model_coneq,output_respsurf]=getRespSurfModel...
        (x_RBF_list,fval_RBF_list,con_RBF_list,coneq_RBF_list);
    object_function_surrogate=output_respsurf.object_function_surrogate;
    nonlcon_function_surrogate=output_respsurf.nonlcon_function_surrogate;
    
    % get local infill point
    % obtian total constraint function
    if ~isempty(nonlcon_function_surrogate) || ~isempty(cheapcon_function)
        constraint_function=@(x) totalconFunction...
            (x,nonlcon_function_surrogate,cheapcon_function,nonlcon_torlance);
    else
        constraint_function=[];
    end
    fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp');
    x_local_infill=fmincon(object_function_surrogate,x_RBF_center,[],[],[],[],...
        low_bou,up_bou,constraint_function,fmincon_options);
    
    [x_local_infill,fval_local_infill,con_local_infill,coneq_local_infill,NFE_p]=dataLibraryUpdataProtect...
        (data_library_name,model_function,x_local_infill',...
        x_list,fval_list,con_list,coneq_list,...
        low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
    
    x_list=[x_list;x_local_infill];
    fval_list=[fval_list;fval_local_infill];
    if ~isempty(con_list)
        con_list=[con_list;con_local_infill];
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_local_infill];
    end
    
    fval_nomlz_list=[fval_nomlz_list;fval_global_infill/fval_max];
    if ~isempty(con_nomlz_list)
        con_nomlz_list=[con_nomlz_list;con_local_infill./con_max_list];
    end
    if ~isempty(con_nomlz_list)
        coneq_nomlz_list=[coneq_nomlz_list;coneq_local_infill./coneq_max_list];
    end
    
    if DRAW_FIGURE_FLAG && variable_number < 3
        interpolationVisualize(respsurf_model_fval,low_bou,up_bou);
        line(x_local_infill(1),x_local_infill(2),fval_local_infill,'Marker','o','color','r')
    end
    
    % find best result to record
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);
    
    if INFORMATION_FLAG
        fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
        fprintf('global x:          %s\n',num2str(x_global_infill));
        fprintf('global value:      %f\n',fval_global_infill);
        fprintf('global violation:  %s  %s\n',num2str(con_global_infill),num2str(coneq_global_infill));
        fprintf('local  x:          %s\n',num2str(x_local_infill));
        fprintf('local  value:      %f\n',fval_local_infill);
        fprintf('local  violation:  %s  %s\n',num2str(con_local_infill),num2str(coneq_local_infill));
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

% save(['result.mat']);
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
    function fval=meritFunction...
            (x,object_function,nonlcon_function,lambda_con,lambda_coneq,miu)
        % penalty function
        % augmented lagrange multiplier method was used
        %
        fval=object_function(x);
        if ~isempty(nonlcon_function)
            [con__,coneq__]=nonlcon_function(x);
            if ~isempty(con__)
                psi=max(con__,-lambda_con/2/miu);
                fval=fval+sum(lambda_con.*psi+miu*psi.*psi);
            end
            if ~isempty(coneq__)
                fval=fval+sum(lambda_coneq.*coneq__+miu*coneq__.*coneq__);
            end
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
        for x_index__=size(x_add_list,1)
            x_updata__=x_add_list(x_index__,:);
            
            % check x_potential if exist in data library
            % if not, updata data libraray
            distance__=sum(((x_updata__-x_list)./(low_bou'-up_bou')).^2,2);
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

function [x_best,fval_best,exitflag,output]=findMinMSP...
    (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function_surrogate,...
    cheapcon_function,nonlcon_torlance)
% find min fval use MSP guideline
% MSP: object_funtion is object_function (generate by surrogate model)
% nonlcon_function generate by surrogate model
% use ga as optimal method
%
if nargin < 6
    cheapcon_function=[];
    if nargin < 5
        nonlcon_function_surrogate=[];
        if nargin < 4
            up_bou=[];
            if nargin < 3
                low_bou=[];
            end
        end
    end
end
% % object function convert to penalty function if input nonlcon function
% if ~isempty(nonlcon_function)
%     object_function=@(x) penaltyFunction(object_function,x,nonlcon_function);
%     constraint_function=cheapcon_function;
% end

% obtian total constraint function
if ~isempty(nonlcon_function_surrogate) || ~isempty(cheapcon_function)
    constraint_function=@(x) totalconFunction...
        (x,nonlcon_function_surrogate,cheapcon_function,nonlcon_torlance);
else
    constraint_function=[];
end

% generate initial population for ga
population_matrix=zeros(max(10,2*variable_number),variable_number);
for population_index=1:size(population_matrix,1)
    x=rand(variable_number,1).*(up_bou-low_bou)+low_bou;
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x);
        while sum([~(con < 0);abs(coneq) < 0])
            x=rand(variable_number,1).*(up_bou-low_bou)+low_bou;
            [con,coneq]=cheapcon_function(x);
        end
    end
    population_matrix(population_index,:)=x';
end

% optiaml
ga_option=optimoptions('ga','FunctionTolerance',1e-2,'ConstraintTolerance',1e-2,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',100,'InitialPopulationMatrix',population_matrix,...
    'display','none');
[x_best,fval_best,exitflag,output]=ga...
    (object_function_surrogate,variable_number,[],[],[],[],low_bou',up_bou',constraint_function,ga_option);
x_best=x_best';
fmincon_option=optimoptions('fmincon','FunctionTolerance',1e-6,'ConstraintTolerance',1e-6,...
    'algorithm','sqp',....
    'display','none');
[x_best,fval_best,exitflag,output]=fmincon...
    (object_function_surrogate,x_best,[],[],[],[],low_bou,up_bou,constraint_function,fmincon_option);

    function fval=penaltyFunction(object_function,x,nonlcon_function)
        fval=object_function(x);
        [con__,coneq__]=nonlcon_function(x);
        fval=fval+10*sum(max(con__,0))+10*sum(abs(coneq__));
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
end
function [x_best,fval_best,LCB_function]=findMinLCB...
    (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function_surrogate,cheapcon_function,...
    object_function_variance,w)
% find min fval use LCB guideline
% LCB: object_funtion is ...
% object_function (generate by surrogate model)...
% subtract sqrt(object_var_function(x)) (generate by surrogate model)
% nonlcon_function generate by surrogate model
% use ga as optimal method
%
if ~isempty(nonlcon_function_surrogate) || ~isempty(cheapcon_function)
    totalcon_function=@(x) totalconFunction(x,nonlcon_function_surrogate,cheapcon_function);
else
    totalcon_function=[];
end
if nargin < 12
    w=1;
end

% generate initial population for ga
population_matrix=zeros(2*variable_number,variable_number);
for population_index=1:2*variable_number
    x=rand(variable_number,1).*(up_bou-low_bou)+low_bou;
    if ~isempty(totalcon_function)
        [con,coneq]=totalcon_function(x);
        while sum(~(con < 0))
            x=rand(variable_number,1).*(up_bou-low_bou)+low_bou;
            [con,coneq]=totalcon_function(x);
        end
    end
    population_matrix(population_index,:)=x';
end

% optiaml
ga_option=optimoptions('ga','FunctionTolerance', 1e-2,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',50,'InitialPopulationMatrix',population_matrix,...
    'display','none');

LCB_function=@(x) object_function_surrogate(x)-w*sqrt(object_function_variance(x));
[x_best,fval_best,~,~]=ga...
    (LCB_function,variable_number,[],[],[],[],low_bou',up_bou',totalcon_function,ga_option);

x_best=x_best';
    function [con,coneq]=totalconFunction(x,nonlcon_function,cheapcon_function)
        con=[];
        coneq=[];
        if ~isempty(nonlcon_function)
            [expencon,expenconeq]=nonlcon_function(x);
            con=[con;expencon];
            coneq=[coneq;expenconeq];
        end
        if ~isempty(cheapcon_function)
            [expencon,expenconeq]=cheapcon_function(x);
            con=[con;expencon];
            coneq=[coneq;expenconeq];
        end
    end
end
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

function [kriging_model_fval,kriging_model_con,kriging_model_coneq,output]=getKrigingModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create kriging model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%

kriging_model_fval=interpolationKrigingPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    kriging_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'base_function_list',[],'fval_regression',[],'covariance',[],'inv_convariance',[],...
        'theta',[],'beta',[],'gama',[],'sigma_sq',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'predict_function',[]);
    kriging_model_con=repmat(kriging_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        kriging_model_con(con_index)=interpolationKrigingPreModel...
            (x_list,con_list(:,con_index));
    end
else
    kriging_model_con=[];
end

if ~isempty(coneq_list)
    kriging_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'base_function_list',[],'fval_regression',[],'covariance',[],'inv_convariance',[],...
        'theta',[],'beta',[],'gama',[],'sigma_sq',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'predict_function',[]);
    kriging_model_coneq=repmat(kriging_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        kriging_model_coneq(coneq_index)=interpolationKrigingPreModel...
            (x_list,coneq_list(:,coneq_index),theta_coneq,base_function_list);
    end
else
    kriging_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,kriging_model_fval);
object_function_variance=@(predict_x) objectFunctionVariance(predict_x,kriging_model_fval);
if isempty(kriging_model_con) && isempty(kriging_model_coneq)
    nonlcon_function_surrogate=[];
    nonlcon_function_variance=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,kriging_model_con,kriging_model_coneq);
    nonlcon_function_variance=@(predict_x) nonlconFunctionVariance(predict_x,kriging_model_con,kriging_model_coneq);
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
            (predict_x,kriging_model_fval)
        [fval,~]=kriging_model_fval.predict_function(predict_x);
    end
    function fval_var=objectFunctionVariance...
            (predict_x,kriging_model_fval)
        [~,fval_var]=kriging_model_fval.predict_function(predict_x);
    end
    function [con,coneq]=nonlconFunctionSurrogate...
            (predict_x,kriging_model_con,kriging_model_coneq)
        if isempty(kriging_model_con)
            con=[];
        else
            con=zeros(length(kriging_model_con),1);
            for con_index__=1:length(kriging_model_con)
                [con(con_index__),~]=kriging_model_con(con_index__).predict_function....
                    (predict_x);
            end
        end
        if isempty(kriging_model_coneq)
            coneq=[];
        else
            coneq=zeros(length(kriging_model_coneq),1);
            for coneq_index__=1:length(kriging_model_con)
                [coneq(coneq_index__),~]=kriging_model_coneq(coneq_index__).predict_function...
                    (predict_x);
            end
        end
    end
    function [con_var,coneq_var]=nonlconFunctionVariance...
            (predict_x,kriging_model_con,kriging_model_coneq)
        if isempty(kriging_model_con)
            con_var=[];
        else
            con_var=zeros(length(kriging_model_con),1);
            for con_index__=1:length(kriging_model_con)
                [~,con_var(con_index__)]=kriging_model_con(con_index__).predict_function....
                    (predict_x);
            end
        end
        if isempty(kriging_model_coneq)
            coneq_var=[];
        else
            coneq_var=zeros(length(kriging_model_coneq),1);
            for coneq_index__=1:length(kriging_model_con)
                [~,coneq_var(coneq_index__)]=kriging_model_coneq(coneq_index__).predict_function...
                    (predict_x);
            end
        end
    end
end
function kriging_model=interpolationKrigingPreModel...
    (X,Y,theta,base_function_list)
% prepare model, optimal theta, load result and calculation parameter
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a kriging model, include X, Y, base_function_list
% and predict_function
% theta beta gama sigma_sq is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
if nargin < 4
    base_functon=@(x) 1;
    base_function_list={base_functon};
    if nargin < 3
        theta=[];
        if nargin < 2
            error('interpolationKrigingPreModel: lack input');
        end
    end
end
[x_number,variable_number]=size(X);
if isempty(theta)
    theta=ones(variable_number,1);
end

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__ = find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__) = 1; end
index__ = find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__) = 1; end
X_nomlz = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_nomlz = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

% initial X_dis_sq
X_dis_sq=zeros(x_number,x_number,variable_number);
for rank=1:x_number
    for colume=1:rank-1
        X_dis_sq(rank,colume,:)=X_dis_sq(colume,rank,:);
    end
    for colume=rank:x_number
        X_dis_sq(rank,colume,:)=(X_nomlz(rank,:)-X_nomlz(colume,:)).^2;
    end
end

% optimal to get hyperparameter
options = optimoptions(@fmincon,'Display','off');
low_bou_kriging=1e-1*ones(variable_number,1);
up_bou_kriging=20*ones(variable_number,1);
object_function_kriging=@(theta) objectFunctionKriging...
    (X_dis_sq,X_nomlz,Y_nomlz,x_number,variable_number,theta,base_function_list);
theta=fmincon...
    (object_function_kriging,theta,[],[],[],[],low_bou_kriging,up_bou_kriging,[],options);

% get parameter
[covariance,inv_convariance,fval_reg,beta,sigma_sq]=interpolationKriging...
    (X_dis_sq,X_nomlz,Y_nomlz,x_number,variable_number,theta,base_function_list);
gama=inv_convariance*(Y_nomlz-fval_reg*beta);

% initialization predict function
predict_function=@(predict_x) interpolationKrigingPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    theta,beta,gama,sigma_sq,...
    inv_convariance,fval_reg,base_function_list,predict_x);

kriging_model.X=X;
kriging_model.Y=Y;
kriging_model.X_normalize=X_nomlz;
kriging_model.Y_normalize=Y_nomlz;
kriging_model.base_function_list=base_function_list;
kriging_model.fval_regression=fval_reg;
kriging_model.covariance=covariance;
kriging_model.inv_convariance=inv_convariance;

kriging_model.theta=theta;
kriging_model.beta=beta;
kriging_model.gama=gama;
kriging_model.sigma_sq=sigma_sq;
kriging_model.aver_X=aver_X;
kriging_model.stdD_X=stdD_X;
kriging_model.aver_Y=aver_Y;
kriging_model.stdD_Y=stdD_Y;

kriging_model.predict_function=predict_function;

    function sigma_sq=objectFunctionKriging...
            (X_dis_sq,X,Y,x_number,variable_number,theta,base_function_list)
        % function to minimize sigma_sq
        %
        [~,~,~,~,sigma_sq]=interpolationKriging...
            (X_dis_sq,X,Y,x_number,variable_number,theta,base_function_list);
    end
    function [covariance,inv_convariance,fval_reg,beta,sigma_sq]=interpolationKriging...
            (X_dis_sq,X,Y,x_number,variable_number,theta,base_function_list)
        % total riging interpolation function
        % input X, Y as initial data, theta and base function
        % output beta and gama, xita
        %
        % Copyright 2022 Adel
        %
        % get covariance matrix
        covariance=zeros(x_number,x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                covariance(rank_index__,colume_index__)=...
                    covariance(colume_index__,rank_index__);
            end
            covariance(rank_index__,rank_index__)=1;
            for colume_index__=rank_index__+1:x_number
                temp__=0;
                for variable_index=1:variable_number
                    temp__=temp__+...
                        X_dis_sq(rank_index__,colume_index__,variable_index)*...
                        theta(variable_index);
                end
                covariance(rank_index__,colume_index__)=...
                    exp(-temp__);
            end
        end
        base_function_number__=size(base_function_list,1);
        
        % stabilize
        covariance=covariance+eye(x_number)*1e-6;
        
        % F is base funcion fval of data_point_x
        fval_reg=zeros(x_number,base_function_number__);
        for base_function_index__=1:base_function_number__
            base_function=base_function_list{base_function_index__};
            for x_index__=1:x_number
                fval_reg(x_index__,base_function_index__)=...
                    base_function(X(x_index__,:));
            end
        end
        
        % coefficient calculation
        inv_convariance=inv(covariance);
        beta=(fval_reg'*inv_convariance*fval_reg)\fval_reg'*inv_convariance*Y;
        sigma_sq=(Y-fval_reg*beta)'*inv_convariance*(Y-fval_reg*beta)/x_number;
    end
    function [predict_fval,predict_variance]=interpolationKrigingPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            theta,beta,gama,sigma_sq,...
            inv_convariance,fval_reg,base_function_list,predict_x)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        % Copyright 2022 Adel
        %
        if size(predict_x,1) > 1
            predict_x=predict_x';
        end
        
        % normalize data
        predict_x=(predict_x-aver_X)./stdD_X;
        
        % predict value
        base_function_number__=size(base_function_list,1);
        predict_cov__=exp(-(predict_x-X_nomlz).^2*theta);
        predict_fval_reg__=zeros(1,base_function_number__);
        for base_function_index=1:base_function_number__
            base_function=base_function_list{base_function_index};
            predict_fval_reg__(1,base_function_index)=base_function(predict_x');
        end
        predict_fval=predict_fval_reg__'*beta+predict_cov__'*gama;
        
        % predict variance
        u__=fval_reg'*inv_convariance*predict_cov__-predict_fval_reg__;
        predict_variance=sigma_sq*...
            (1+u__'/(fval_reg'*inv_convariance*fval_reg)*u__+...
            -predict_cov__'*inv_convariance*predict_cov__);
        
        % normalize data
        predict_fval=predict_fval*stdD_Y+aver_Y;
        predict_variance=predict_variance*stdD_Y*stdD_Y;
    end
end

function [ERBF_model_fval,ERBF_model_con,ERBF_model_coneq,output]=getEnsemleRadialBasisModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create kriging model and function
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
%
ERBF_model_fval=interpolationEnsemleRadialBasisPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    ERBF_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'basis_function_list',[],'c_list',[],'beta_list',[],...
        'rdibas_matrix_list',[],'inv_rdibas_matrix_list',[],'model_error_list',[],'w',[],...
        'predict_function',[]);
    ERBF_model_con=repmat(ERBF_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        ERBF_model_con(con_index)=interpolationEnsemleRadialBasisPreModel...
            (x_list,con_list(:,con_index));
    end
else
    ERBF_model_con=[];
end

if ~isempty(coneq_list)
    ERBF_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'basis_function_list',[],'c_list',[],'beta_list',[],...
        'rdibas_matrix_list',[],'inv_rdibas_matrix_list',[],'model_error_list',[],'w',[],...
        'predict_function',[]);
    ERBF_model_coneq=repmat(ERBF_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        ERBF_model_coneq(coneq_index)=interpolationEnsemleRadialBasisPreModel...
            (x_list,coneq_list(:,coneq_index));
    end
else
    ERBF_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,ERBF_model_fval);
if isempty(ERBF_model_con) && isempty(ERBF_model_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,ERBF_model_con,ERBF_model_coneq);
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
function ensemleradialbasis_model=interpolationEnsemleRadialBasisPreModel...
    (X,Y)
% get ensemle radial basis function interpolation model function version 1
% use cubic interpolation optimal to decrese time use
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
% beta is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
[x_number,variable_number]=size(X);

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index = find(stdD_X == 0);
if ~isempty(index),stdD_X(index)=1;end
index = find(stdD_Y == 0);
if ~isempty(index),stdD_Y(index)=1;end
X_nomlz = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_nomlz = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

c_initial=(prod(max(X_nomlz)-min(Y_nomlz))/x_number)^(1/variable_number);

% initial sq of x and x
x_sq_matrix=zeros(x_number);
for rank_index=1:x_number
    for colume_index=1:rank_index-1
        x_sq_matrix(rank_index,colume_index)=...
            x_sq_matrix(colume_index,rank_index);
    end
    for colume_index=rank_index:x_number
        x_sq_matrix(rank_index,colume_index)=...
            sum((X_nomlz(rank_index,:)'-X_nomlz(colume_index,:)').^2);
    end
end

option_fmincon=optimoptions('fmincon','display','none','TolFun',1e-6,...
    'SpecifyObjectiveGradient',true);

% linear kernal function
basis_function_linear=@(x_sq,c) sqrt(x_sq)+c;
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) ones(x_number);
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_linear,c,rdibas_matrix_gradient_function);
[linear_c_backward,fval_backward,NFE]=optimalCubicInterpolation...
    (object_function,-1e2,-1e2,1e2,1e-3);
[linear_c_forward,fval_forward,NFE]=optimalCubicInterpolation...
    (object_function,1e2,-1e2,1e2,1e-3);
if fval_forward < fval_backward
    c_linear=linear_c_forward;
else
    c_linear=linear_c_backward;
end

% gauss kernal function
basis_function_gauss=@(x_sq,c) exp(-c*x_sq);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) -x_sq_matrix.*rdibas_matrix;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_gauss,c,rdibas_matrix_gradient_function);
[c_gauss,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% spline kernal function
basis_function_spline=@(x_sq,c) x_sq*log(x_sq*c+1e-3);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) x_sq_matrix.*x_sq_matrix./(x_sq_matrix*c+1e-3);
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_spline,c,rdibas_matrix_gradient_function);
[c_spline,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% triple kernal function
basis_function_triple=@(x_sq,c) (sqrt(x_sq)+c)^3;
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) 3*(sqrt(x_sq_matrix)+c).^3;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_triple,c,rdibas_matrix_gradient_function);
[c_triple,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% multiquadric kernal function
basis_function_multiquadric=@(x_sq,c) sqrt(x_sq+c);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) 0.5./rdibas_matrix;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_multiquadric,c,rdibas_matrix_gradient_function);
[c_binomial,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% inverse multiquadric kernal function
basis_function_inverse_multiquadric=@(x_sq,c) 1/sqrt(x_sq+c);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) -0.5*rdibas_matrix.^3;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_inverse_multiquadric,c,rdibas_matrix_gradient_function);
[c_inverse_binomial,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);
% c_initial=1;
% [fval,gradient]=object_function(c_initial)
% [fval_diff,gradient_diff]=differ(object_function,c_initial)
% drawFunction(object_function,1e-1,10);

% generate total model
basis_function_list={
    basis_function_linear;
    basis_function_gauss;
    basis_function_spline;
    basis_function_triple;
    basis_function_multiquadric;
    basis_function_inverse_multiquadric;};
c_list=[
    c_linear;
    c_gauss;
    c_spline;
    c_triple;
    c_binomial;
    c_inverse_binomial;
    ];

model_number=size(basis_function_list,1);
beta_list=zeros(x_number,1,model_number);
rdibas_matrix_list=zeros(x_number,x_number,model_number);
inv_rdibas_matrix_list=zeros(x_number,x_number,model_number);

% calculate model matrix and error
model_error_list=zeros(model_number,x_number);
for model_index=1:model_number
    basis_function=basis_function_list{model_index};
    c=c_list(model_index);
    [beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
        (x_sq_matrix,Y_nomlz,x_number,basis_function,c);
    beta_list(:,:,model_index)=beta;
    rdibas_matrix_list(:,:,model_index)=rdibas_matrix;
    inv_rdibas_matrix_list(:,:,model_index)=inv_rdibas_matrix;
    
    model_error_list(model_index,:)=(beta./...
        diag(inv_rdibas_matrix))';
end

% calculate weight of each model
C=model_error_list*model_error_list';
eta=trace(C)/x_number;
w=(C+eta*eye(model_number))\ones(model_number,1)/...
    (ones(1,model_number)/(C+eta*eye(model_number))*ones(model_number,1));
while min(w) < -0.05
    eta=eta*10;
    w=(C+eta*eye(model_number))\ones(model_number,1)/...
        (ones(1,model_number)/(C+eta*eye(model_number))*ones(model_number,1));
end

% initialization predict function
predict_function=@(predict_x) interpolationEnsemleRadialBasisPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    model_number,beta_list,basis_function_list,c_list,...
    w,predict_x);

ensemleradialbasis_model.X=X;
ensemleradialbasis_model.Y=Y;
ensemleradialbasis_model.X_normalize=X_nomlz;
ensemleradialbasis_model.Y_normalize=Y_nomlz;
ensemleradialbasis_model.aver_X=aver_X;
ensemleradialbasis_model.stdD_X=stdD_X;
ensemleradialbasis_model.aver_Y=aver_Y;
ensemleradialbasis_model.stdD_Y=stdD_Y;

ensemleradialbasis_model.basis_function_list=basis_function_list;
ensemleradialbasis_model.c_list=c_list;
ensemleradialbasis_model.beta_list=beta_list;
ensemleradialbasis_model.rdibas_matrix_list=rdibas_matrix_list;
ensemleradialbasis_model.inv_rdibas_matrix_list=inv_rdibas_matrix_list;
ensemleradialbasis_model.model_error_list=model_error_list;
ensemleradialbasis_model.w=w;

ensemleradialbasis_model.predict_function=predict_function;

    function [fval,gradient]=objectFunctionRadiabasis....
            (x_sq_matrix,Y,x_number,basis_function,c,rdibas_matrix_gradient_function)
        % MSE_CV function, simple approximation to RMS
        % basis_function input is c and x_sq
        %
        [beta__,rdibas_matrix__,inv_rdibas_matrix__]=interpolationRadialBasis...
            (x_sq_matrix,Y,x_number,basis_function,c);
        fval=0;
        U=zeros(x_number,1);
        for x_index=1:x_number
            U(x_index)=...
                beta__(x_index)/inv_rdibas_matrix__(x_index,x_index);
            fval=fval+(U(x_index))^2;
        end
        
        % calculate gradient
        if nargout > 1
            inv_rdibas_matrix_gradient=-inv_rdibas_matrix__*...
                rdibas_matrix_gradient_function...
                (x_number,x_sq_matrix,rdibas_matrix__,c)*inv_rdibas_matrix__;
            U_gradient=zeros(x_number,1);
            I=eye(x_number);
            for x_index=1:x_number
                U_gradient(x_index)=(I(x_index,:)*inv_rdibas_matrix_gradient*Y)/...
                    inv_rdibas_matrix__(x_index,x_index)-...
                    beta__(x_index)*(I(x_index,:)*inv_rdibas_matrix_gradient*I(:,x_index))/...
                    inv_rdibas_matrix__(x_index,x_index)^2;
            end
            
            gradient=2*sum(U.*U_gradient);
        end
    end
    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
            (x_sq_matrix,Y,x_number,basis_function,c)
        % interpolation polynomial responed surface core function
        % calculation beta
        rdibas_matrix=zeros(x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                rdibas_matrix(rank_index__,colume_index__)=...
                    rdibas_matrix(colume_index__,rank_index__);
            end
            for colume_index__=rank_index__:x_number
                rdibas_matrix(rank_index__,colume_index__)=...
                    basis_function(x_sq_matrix(rank_index__,colume_index__),c);
            end
        end
        
        % stabilize matrix
        rdibas_matrix=rdibas_matrix+eye(x_number)*1e-6;
        
        if rcond(rdibas_matrix) < eps || isnan(rdibas_matrix(1))
            disp('error');
        end
        
        inv_rdibas_matrix=inv(rdibas_matrix);
        beta=inv_rdibas_matrix*Y;
    end
    function [x_best,favl_best,NFE,output]=optimalCubicInterpolation...
            (object_function,x_initial,low_bou,up_bou,torlance,iteration_max)
        % cubic interpolation optimization, should provide fval and gradient
        % only work for one best(convex)
        %
        if nargin < 6
            iteration_max=[];
            if nargin < 5
                torlance=[];
                if nargin < 4
                    up_bou=[];
                    if nargin < 3
                        low_bou=[];
                        if nargin < 2
                            error('lack x initial');
                        end
                    end
                end
            end
        end
        
        INFORMATION_FLAG=0; % whether show information
        
        draw_range=0.001;
        draw_interval=draw_range*0.02;
        DRAW_FLAG=0;
        
        if isempty(iteration_max)
            iteration_max=10*length(x_initial);
        end
        
        if isempty(torlance)
            torlance=1e-6;
        end
        
        x=x_initial;
        done=0;
        iteration=0;
        NFE=0;
        result_x_list=[];
        result_fval_list=[];
        
        % decide which turn to search
        [fval,gradient]=object_function(x);NFE=NFE+1;
        result_x_list=[result_x_list;x];
        result_fval_list=[result_fval_list;fval];
        if gradient < -torlance
            direction=1;
        elseif gradient > torlance
            direction=-1;
        else
            done=1;
            x_best=x;
            favl_best=fval;
        end
        
        x_old=x;
        fval_old=fval;
        gradient_old=gradient;
        iteration=iteration+1;
        
        % move forward to first point
        if ~done
            x=x_old+direction*0.01;
            if x > up_bou
                x=up_bou;
            elseif x < low_bou
                x=low_bou;
            end
            [fval,gradient]=object_function(x);NFE=NFE+1;
            result_x_list=[result_x_list;x];
            result_fval_list=[result_fval_list;fval];
            quit_flag=judgeQuit...
                (x,x_old,fval,fval_old,gradient,torlance,iteration,iteration_max);
            if quit_flag
                done=1;
                x_best=x;
                favl_best=fval;
            end
            iteration=iteration+1;
        end
        
        % main loop for cubic interpolation
        while ~done
            
            x_base=x_old;
            x_relative=x/x_old;
            interp_matrix=[1,1,1,1;
                3,2,1,0;
                x_relative^3,x_relative^2,x_relative,1;
                3*x_relative^2,2*x_relative,1,0];
            
            if rcond(interp_matrix) < eps
                disp('error');
            end
            
            interp_value=[fval_old;gradient_old*x_base;fval;gradient*x_base];
            [x_inter_rel,coefficient_cubic]=minCubicInterpolate(interp_matrix,interp_value);
            x_inter=x_inter_rel*x_base;
            
            if DRAW_FLAG
                x_draw=1:direction*draw_interval:direction*draw_range;
                x_draw=x_draw/x_base;
                line(x_draw*x_base,coefficient_cubic(1)*x_draw.^3+coefficient_cubic(2)*x_draw.^2+...
                    coefficient_cubic(3)*x_draw+coefficient_cubic(4));
            end
            
            % limit search space, process constraints
            if x_inter > up_bou
                x_inter=up_bou;
            elseif x_inter < low_bou
                x_inter=low_bou;
            end
            
            [fval_inter,gradient_inter]=object_function(x_inter);NFE=NFE+1;
            
            % only work for one best(convex)
            % three situation discuss
            if gradient < 0
                x_old=x;
                fval_old=fval;
                gradient_old=gradient;
            else
                if gradient_inter < 0
                    x_old=x;
                    fval_old=fval;
                    gradient_old=gradient;
                end
            end
            
            x=x_inter;
            fval=fval_inter;
            gradient=gradient_inter;
            
            quit_flag=judgeQuit...
                (x,x_old,fval,fval_old,gradient,torlance,iteration,iteration_max);
            if quit_flag
                done=1;
                x_best=x;
                favl_best=fval;
            end
            
            result_x_list=[result_x_list;x];
            result_fval_list=[result_fval_list;fval];
            iteration=iteration+1;
        end
        output.result_x_list=result_x_list;
        output.result_fval_list=result_fval_list;
        
        function [lamada,coefficient_cubic]=minCubicInterpolate(interpolate_matrix,interpolate_value)
            % calculate min cubic curve
            %
            coefficient_cubic=interpolate_matrix\interpolate_value;
            
            temp_sqrt=4*coefficient_cubic(2)^2-12*coefficient_cubic(1)*coefficient_cubic(3);
            if temp_sqrt>=0
                temp_lamada=-coefficient_cubic(2)/3/coefficient_cubic(1)+...
                    sqrt(temp_sqrt)/6/coefficient_cubic(1);
                if (temp_lamada*6*coefficient_cubic(1)+2*coefficient_cubic(2))>0
                    lamada=temp_lamada;
                else
                    lamada=-coefficient_cubic(2)/3/coefficient_cubic(1)-...
                        sqrt(temp_sqrt)...
                        /6/coefficient_cubic(1);
                end
            else
                lamada=-coefficient_cubic(2)/3/coefficient_cubic(1);
            end
        end
        function quit_flag=judgeQuit...
                (x,x_old,fval,fval_old,gradient,torlance,iteration,iteration_max)
            quit_flag=0;
            if abs(fval-fval_old)/fval_old < torlance
                quit_flag=1;
            end
            if abs(gradient) < torlance
                quit_flag=1;
            end
            if abs(x-x_old) < 1e-4
                quit_flag=1;
            end
            if iteration >= iteration_max
                quit_flag=1;
            end
        end
    end
    function predict_y=interpolationEnsemleRadialBasisPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            model_number,beta_list,basis_function_list,c_list,...
            w,predict_x)
        % ensemle radial basis function interpolation predict function
        % input predict_x and radialbasis_model model
        if size(predict_x,2) > 1
            predict_x=predict_x';
        end
        
        % normalize data
        predict_x_nomlz=(predict_x-aver_X')./stdD_X';
        
        % calculate each sub model predict fval and get predict_y
        predict_y_nomlz=0;
        for model_index__=1:model_number
            basis_function__=basis_function_list{model_index__};
            c__=c_list(model_index__);
            beta__=beta_list(:,:,model_index__);
            predict_y_nomlz_sub=interpolationRadialBasisPredictor...
                (X_nomlz,beta__,basis_function__,c__,predict_x_nomlz);
            predict_y_nomlz=predict_y_nomlz+w(model_index__)*predict_y_nomlz_sub;
        end
        
        % normalize data
        predict_y=predict_y_nomlz*stdD_Y+aver_Y;
        
        function predict_y_nomlz=interpolationRadialBasisPredictor...
                (X_nomlz,beta,basis_function,c,predict_x_nomlz)
            % radial basis function interpolation predict function
            % input predict_x and radialbasis_model model
            % predict_x is row vector
            % output the predict value
            [x_number__,~]=size(X_nomlz);
            
            % predict value
            X_inner_product__=zeros(x_number__,1);
            for x_index__=1:x_number__
                X_inner_product__(x_index__)=...
                    basis_function(sum((X_nomlz(x_index__,:)'-predict_x_nomlz).^2),c);
            end
            
            % predict variance
            predict_y_nomlz=beta'*X_inner_product__;
        end
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
radialbasis_model_fval=interpolationRadialBasisPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    radialbasis_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'radialbasis_matrix',[],'inv_radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_con=repmat(radialbasis_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        radialbasis_model_con(con_index)=interpolationRadialBasisPreModel...
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
        radialbasis_model_coneq(coneq_index)=interpolationRadialBasisPreModel...
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
function radialbasis_model=interpolationRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interpolation pre model function version 1
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

[beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
    (X_nomlz,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function=@(predict_x) interpolationRadialBasisPredictor...
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

    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
            (X,Y,basis_function,x_number)
        % interpolation polynomial responed surface core function
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
    function [predict_y]=interpolationRadialBasisPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            beta,basis_function,predict_x)
        % radial basis function interpolation predict function
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

function [respsurf_model_fval,respsurf_model_con,respsurf_model_coneq,output]=getRespSurfModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create respsurf model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
respsurf_model_fval=interpolationRespSurfPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    respsurf_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'predict_function',[]);
    respsurf_model_con=repmat(respsurf_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        respsurf_model_con(con_index)=interpolationRespSurfPreModel...
            (x_list,con_list(:,con_index));
    end
else
    respsurf_model_con=[];
end

if ~isempty(coneq_list)
    respsurf_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'predict_function',[]);
    respsurf_model_coneq=repmat(respsurf_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        respsurf_model_coneq(coneq_index)=interpolationRespSurfPreModel...
            (x_list,coneq_list(:,coneq_index));
    end
else
    respsurf_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,respsurf_model_fval);
if isempty(respsurf_model_con) && isempty(respsurf_model_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,respsurf_model_con,respsurf_model_coneq);
end

output.object_function_surrogate=object_function_surrogate;
output.nonlcon_function_surrogate=nonlcon_function_surrogate;
output.x_list=x_list;
output.fval_list=fval_list;
output.con_list=con_list;
output.coneq_list=coneq_list;

    function fval=objectFunctionSurrogate...
            (predict_x,respsurf_model_fval)
        fval=respsurf_model_fval.predict_function(predict_x);
    end
    function [con,coneq]=nonlconFunctionSurrogate...
            (predict_x,respsurf_model_con,respsurf_model_coneq)
        if isempty(respsurf_model_con)
            con=[];
        else
            con=zeros(length(respsurf_model_con),1);
            for con_index__=1:length(respsurf_model_con)
                con(con_index__)=respsurf_model_con...
                    (con_index__).predict_function(predict_x);
            end
        end
        if isempty(respsurf_model_coneq)
            coneq=[];
        else
            coneq=zeros(length(respsurf_model_coneq),1);
            for coneq_index__=1:length(respsurf_model_con)
                coneq(coneq_index__)=respsurf_model_coneq...
                    (coneq_index__).predict_function(predict_x);
            end
        end
    end
end
function respsurf_model=interpolationRespSurfPreModel(X,Y)
% polynomial response surface interpolation pre model function
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
% beta is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
x_number=size(X,1);
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

beta=interpolationRespSurf(X_nomlz,Y_nomlz);

% initialization predict function
predict_function=@(predict_x) interpolationRespSurfPredictor...
    (aver_X,stdD_X,aver_Y,stdD_Y,beta,predict_x);

respsurf_model.X=X;
respsurf_model.Y=Y;
respsurf_model.X_normalize=X_nomlz;
respsurf_model.Y_normalize=Y_nomlz;
respsurf_model.beta=beta;
respsurf_model.aver_X=aver_X;
respsurf_model.stdD_X=stdD_X;
respsurf_model.aver_Y=aver_Y;
respsurf_model.stdD_Y=stdD_Y;
respsurf_model.predict_function=predict_function;

    function beta=interpolationRespSurf(X,Y)
        % interpolation polynomial responed surface core function
        % calculation beta
        %
        [x_number__,variable_number__]=size(X);
        
        X_inner__=zeros(x_number__,(variable_number__+1)*(variable_number__+2)/2);
        x_cross__=zeros(1,(variable_number__-1)*variable_number__/2);
        for x_index=1:x_number__
            x__=X(x_index,:);
            x_sqrt__=x__.^2;
            cross_index__=1;
            for i_index=1:variable_number__
                for j_index=i_index+1:variable_number__
                    x_cross__(cross_index__)=x__(i_index)*x__(j_index);
                    cross_index__=cross_index__+1;
                end
            end
            X_inner__(x_index,:)=[1,x__,x_sqrt__,x_cross__];
        end
        
        X_inter_X_inter__=X_inner__'*X_inner__;
        X_inter_X_inter__=X_inter_X_inter__+eye((variable_number__+1)*(variable_number__+2)/2)*1e-3;
        beta=X_inter_X_inter__\X_inner__'*Y;
    end
    function [predict_y]=interpolationRespSurfPredictor...
            (aver_X,stdD_X,aver_Y,stdD_Y,beta,predict_x)
        % polynomial response surface interpolation predict function
        % input predict_x and respsurf_model model
        % predict_x is row vector
        % output the predict value
        %
        % Copyright 2022 Adel
        %
        if size(predict_x,1) > 1
            predict_x=predict_x';
        end
        variable_number__=size(predict_x,2);
        
        % normalize data
        predict_x=(predict_x-aver_X)./stdD_X;
        
        % predict value
        x_cross__=zeros(1,(variable_number__-1)*variable_number__/2);
        x_sqrt__=predict_x.^2;
        cross_index__=1;
        for i_index=1:variable_number__
            for j_index=i_index+1:variable_number__
                x_cross__(cross_index__)=predict_x(i_index)*predict_x(j_index);
                cross_index__=cross_index__+1;
            end
        end
        x_inter=[1,predict_x,x_sqrt__,x_cross__];
        
        % predict variance
        predict_y=x_inter*beta;
        
        % normalize data
        predict_y=predict_y*stdD_Y+aver_Y;
    end
end

function interpolationVisualize...
    (model,low_bou,up_bou,figure_handle)
% visualization polynamial respond surface model
% figrue is 100
%
% Copyright 2022 Adel
%
if nargin < 4
    figure_handle=figure(101);
    if nargin < 3
        up_bou=[];
        if nargin < 2
            low_bou=[];
        end
    end
end
if size(low_bou,1) ~= size(low_bou,1)
    error('interpolationRadialBasisVisualize: boundary incorrect');
end
if size(low_bou,1) > 2
    error('interpolationRadialBasisVisualize: dimension large than two');
end

delete(figure_handle.Children);
axes_handle=axes(figure_handle);

x_list=model.X;
y_list=model.Y;
predict_function=model.predict_function;

% get boundary
if isempty(low_bou)
    low_bou=min(x_list,[],1)';
end
if isempty(up_bou)
    up_bou=max(x_list,[],1)';
end

grid_number=100;
d_bou=(up_bou-low_bou)/grid_number;

if size(x_list,2) == 1
    predict_result=zeros(grid_number+1,1);
    X_draw=low_bou:d_bou:(low_bou+grid_number*d_bou);
    for x_index=1:grid_number+1
        predict_x=(x_index-1).*d_bou+low_bou;
        predict_result(x_index)=predict_function(predict_x);
    end
    line(axes_handle,X_draw,predict_result);
    line(axes_handle,x_list,y_list,'Marker','o','LineStyle','none');
    xlabel('X');
    ylabel('Y');
elseif size(x_list,2) == 2
    predict_result=zeros(grid_number+1);
    [X_draw,Y_draw]=meshgrid(low_bou(1):d_bou(1):(low_bou(1)+grid_number*d_bou(1)),...
        low_bou(2):d_bou(2):(low_bou(2)+grid_number*d_bou(2)));
    for x_index=1:grid_number+1
        for y_index=1:grid_number+1
            predict_x=([x_index,y_index]-1).*d_bou'+low_bou';
            predict_result(y_index,x_index)=predict_function(predict_x);
        end
    end
    surf(axes_handle,X_draw,Y_draw,predict_result,'FaceAlpha',0.5,'EdgeColor','none');
    line(axes_handle,x_list(:,1),x_list(:,2),y_list,'Marker','o','LineStyle','none');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
end

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
            judge=sum(x < low_bou')+sum(x > up_bou');
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
        low_bou=zeros(variable_number,1);
        up_bou=ones(variable_number,1);
    end
end
iteration_max=10*variable_number;

% check x_exist_list if meet boundary
if ~isempty(X_exist)
    index=find(X_exist < low_bou');
    index=[index,find(X_exist > up_bou')];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list range error');
    end
    if size(X_exist,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
    X_exist_nomlz=(X_exist-low_bou')./(up_bou'-low_bou');
else
    X_exist_nomlz=[];
end

% check x_new_number
x_new_number=sample_number-size(X_exist,1);
if x_new_number < 0
    X=X_exist;
    X_new=[];
    distance_min_nomlz=getMinDistance(X_exist_nomlz);
    return;
end

low_bou_nomlz=zeros(variable_number,1);
up_bou_nomlz=ones(variable_number,1);

% initialize X_new by sample with constraint function
% generate sample sequence latin hypercube
% election sequential method is used
% sample number is total point in area
% default low_bou is 0, up_bou is 1, cheapcon_function is []
% low_bou and up_bou is colume vector
% x in x_exist_list, x_list, supply_x_list is row vector
% x_exist_list should meet bou
% get quasi-feasible point
x_initial_number=100*x_new_number;
if ~isempty(cheapcon_function)
    X_supply_quasi_nomlz=[];
    % check if have enough X_supply_nomlz
    while size(X_supply_quasi_nomlz,1) < 100*x_new_number
        X_supply_initial_nomlz=rand(x_initial_number,variable_number);
        x_index=1;
        while x_index <= size(X_supply_initial_nomlz,1)
            x_supply=X_supply_initial_nomlz(x_index,:)'.*(up_bou-low_bou)+low_bou;
            if cheapcon_function(x_supply) > 0
                X_supply_initial_nomlz(x_index,:)=[];
            else
                x_index=x_index+1;
            end
        end
        X_supply_quasi_nomlz=[X_supply_quasi_nomlz;X_supply_initial_nomlz];
    end
else
    X_supply_quasi_nomlz=rand(x_initial_number,variable_number);
end

% iterate and get final x_supply_list
iteration=0;
x_supply_quasi_number=size(X_supply_quasi_nomlz,1);
distance_min_nomlz__=0;
X_new_nomlz=[];
while iteration <= iteration_max
    % random select x_new_number X to X_trial_nomlz
    x_select_index=randperm(x_supply_quasi_number,x_new_number);
    
    % get distance min itertion X_
    distance_min_iteration=getMinDistanceIter...
        (X_supply_quasi_nomlz(x_select_index,:),X_exist_nomlz);
    
    % if distance_min_iteration is large than last time
    if distance_min_iteration > distance_min_nomlz__
        distance_min_nomlz__=distance_min_iteration;
        X_new_nomlz=X_supply_quasi_nomlz(x_select_index,:);
    end
    
    iteration=iteration+1;
end
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);

% x is nomalize, so constraint function should change
if ~isempty(cheapcon_function)
    cheapcon_function=@(x) ...
        max(max(sample_number*cheapcon_function(x.*(up_bou-low_bou)+low_bou)+1,0),[],1);
end

iteration=0;
while iteration < iteration_max
    %     scatter(X_new_nomlz(:,1),X_new_nomlz(:,2));
    %     bou=[low_bou_nomlz,up_bou_nomlz]';
    %     axis(bou(:));axis equal
    %     radius=1/6;
    %     hold on;
    %     rectangle('Position',[-radius+0.5,-radius+0.5,2*radius,2*radius],'Curvature',[1 1])
    %     hold off;
    
    % change each x place by newton methods
    fval_list=zeros(x_new_number,1);
    gradient_list=zeros(variable_number,x_new_number);
    for x_new_index=1:x_new_number
        % get gradient
        [fval_list(x_new_index,1),gradient_list(:,x_new_index),~]=objectFunctionXPlace...
            (X_new_nomlz(x_new_index,:)',[X_new_nomlz(1:x_new_index-1,:);X_new_nomlz(x_new_index+1:end,:);X_exist_nomlz],...
            variable_number,low_bou_nomlz,up_bou_nomlz,cheapcon_function);
    end
    
    % normalize fval
    fval_list=fval_list/max(fval_list);
    for x_new_index=1:x_new_number
        x=X_new_nomlz(x_new_index,:)'-...
            fval_list(x_new_index,1)*gradient_list(:,x_new_index)/norm(gradient_list(:,x_new_index),2)*...
            distance_min_nomlz*(1-iteration/iteration_max);
        for variable_index=1:variable_number
            if x(variable_index,1) > 1
                x(variable_index,1)=1;
            end
            if x(variable_index,1) < 0
                x(variable_index,1)=0;
            end
        end
        X_new_nomlz(x_new_index,:)=x';
    end
    
    iteration=iteration+1;
end
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);
X_new=X_new_nomlz.*(up_bou'-low_bou')+low_bou';
X=[X_new;X_exist];

    function [fval,gradient,hessian]=objectFunctionXPlace...
            (x,X_surplus,variable_number,low_bou,up_bou,cheapcon_function)
        % function describe distance between X and X_supply
        % x is colume vector and X_surplus is matrix which is num-1 x var
        % low_bou_limit__ and up_bou_limit__ is colume vector
        % variable in colume
        %
        [~,variable_number__]=size(X_surplus);
        X_surplus=X_surplus';
        
        sigma__=10;
        boundary__=0.1^variable_number__;
        
        sign__=((x>X_surplus)-0.5)*2;
        
        xi__=-sigma__*(x-X_surplus).*sign__;
        sum_xi__=sum(xi__,1);
        psi__=sigma__*(low_bou-x);
        zeta__=sigma__*(x-up_bou);
        
        exp_sum_xi__=exp(sum_xi__);
        exp_psi__=exp(psi__);
        exp_zeta__=exp(zeta__);
        
        xi_DF=-sigma__*sign__;
        % sum_xi_DF=sum(xi_DF,2);
        psi_DF=-sigma__*ones(variable_number__,1);
        zeta_DF=sigma__*ones(variable_number__,1);
        
        % get fval
        fval=sum(exp_sum_xi__,2)+...
            sum(boundary__*exp_psi__+...
            boundary__*exp_zeta__,1);
        
        % get gradient
        gradient=sum(exp_sum_xi__.*xi_DF,2)+...
            boundary__*exp_psi__.*psi_DF+...
            boundary__*exp_zeta__.*zeta_DF;
        
        % get hessian
        hessian=exp_sum_xi__.*xi_DF*xi_DF'+...
            diag(boundary__*exp_psi__.*psi_DF.*psi_DF+...
            boundary__*exp_zeta__.*zeta_DF.*zeta_DF);
        
        if ~isempty(cheapcon_function)
            fval_con=cheapcon_function(x);
            fval=fval+fval_con;
            [gradient_con,hessian_con]=differ...
                (cheapcon_function,x,fval_con,variable_number);
            gradient=gradient+gradient_con;
            hessian=hessian+hessian_con;
        end
        
        function [gradient,hessian]=differ(differ_function,x,fval,variable_number,step)
            % differ function to get gradient and hessian
            %
            if nargin < 5
                step=1e-6;
            end
            fval__=zeros(variable_number,2); % backward is 1, forward is 2
            gradient=zeros(variable_number,1);
            hessian=zeros(variable_number);
            
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
            
            % hessian
            for variable_index__=1:variable_number
                hessian(variable_index__,variable_index__)=...
                    (fval__(variable_index__,2)-2*fval+fval__(variable_index__,1))/step/step;
                for variable_index_next__=variable_index__+1:variable_number
                    x_for_for=x;
                    x_for_for(variable_index__)=x_for_for(variable_index__)+step;
                    x_for_for(variable_index_next__)=x_for_for(variable_index_next__)+step;
                    
                    hessian(variable_index__,variable_index_next__)=(...
                        differ_function(x_for_for)-...
                        fval__(variable_index__,2)-fval__(variable_index_next__,2)+...
                        fval...
                        )/step/step;
                    hessian(variable_index_next__,variable_index__)=...
                        hessian(variable_index__,variable_index_next__);
                end
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
    function distance_min__=getMinDistanceIter...
            (x_list__,x_exist_list__)
        % get distance min from x_list
        %
        % sort x_supply_list_initial to decrese distance calculate times
        x_list__=sortrows(x_list__,1);
        [sample_number__,variable_number__]=size(x_list__);
        distance_min__=variable_number__;
        for x_index__=1:sample_number__
            x_curr__=x_list__(x_index__,:);
            x_next_index__=x_index__ + 1;
            % only search in min_distance(x_list had been sort)
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
            for x_exist_index=1:size(x_exist_list__,1)
                x_next__=x_exist_list__(x_exist_index,:);
                distance_temp__=sum((x_next__-x_curr__).^2);
                if distance_temp__ < distance_min__
                    distance_min__ = distance_temp__;
                end
            end
        end
        distance_min__=sqrt(distance_min__);
    end
end

function drawFunction(object_function,low_bou,up_bou,...
    grid_number,Y_min,Y_max,figure_handle)
if nargin < 7
    figure_handle=figure(10);
    if nargin < 6
        Y_max=inf;
        if nargin < 5
            Y_min=-inf;
            if nargin < 4
                grid_number=100;
            end
        end
    end
end
axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end
axes_context=axes_handle.Children;
dimension=size(low_bou,1);

switch dimension
    case 1
        d_bou=(up_bou-low_bou)/grid_number;
        X__=low_bou:d_bou:up_bou;
        fval__=zeros(grid_number+1,1);
        for x_index__=1:(grid_number+1)
            fval__(x_index__)=object_function(X__(x_index__));
        end
        line(axes_handle,X__,fval__);
        xlabel('X');
        ylabel('value');
        
    case 2
        d_bou=(up_bou-low_bou)/grid_number;
        [X__,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
        fval__=zeros(grid_number+1);
        for x_index__=1:grid_number+1
            for y_index__=1:grid_number+1
                predict_x=([x_index__;y_index__]-1).*d_bou+low_bou;
                fval__(x_index__,y_index__)=object_function(predict_x);
            end
        end
        fval__(find(fval__ > Y_max))=Y_max;
        fval__(find(fval__ < Y_min))=Y_min;
        axes_context=[axes_context;surface(X__',Y',fval__,'FaceAlpha',0.5,'EdgeColor','none')];
        axes_handle.set('Children',axes_context);
        xlabel('X');
        ylabel('Y');
        zlabel('value');
end
end