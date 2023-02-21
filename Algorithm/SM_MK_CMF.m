clc;
clear;
close all hidden;

benchmark_function=BenchmarkFunction();

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
cheapcon_function=[];

data_library_name='optimal_data_library.txt';
delete(data_library_name);
delete('result_total.txt');
[x_best,fval_best,NFE,output]=optimalSurrogateMKCMF...
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

function [x_best,fval_best,NFE,output]=optimalSurrogateMKCMF...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance)
% surrogate base optimal use radias base function method version 0
% use SVM to get interest point
% FS_FCM to get interest point center point
% and updata interest space
% all function exchange x should be colume vector
% x_list is x_number x variable_number matrix
% both nonlcon_function and cheapcon_function format is [con,coneq]
% model_function should output fval, format is [fval,con,coneq]
% con or coneq can be colume vector if there was more than one constrain
%
% Copyright Adel 2022.10
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

if isempty(iteration_max)
    iteration_max=100;
end

sample_number_initial=min((variable_number+1)*(variable_number+2)/2,5*variable_number);
sample_number_iteration=variable_number;
sample_number_data=10*sample_number_initial;
eta=1/variable_number; % space decrease coefficient

% parameter
scale_SVM=variable_number;
penalty_SVM=100;
% kernal_function_FS_FCM=@(sq) exp(-sq/2/0.1^2);
m=2; % clustering parameter

filter_torlance=1e-3;
protect_range=1e-4;

data_library_name='optimal_data_library';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;

result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);
x_data_list=lhsdesign(sample_number_data,variable_number).*...
    (up_bou'-low_bou')+low_bou';

% if do not input model_function, generate model_function
if nargin < 7 || isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end

iteration=iteration+1;

% step 2
% use latin hypercube method to get initial sample x_list
[~,x_updata_list,~]=getLatinHypercube...
    (sample_number_initial,variable_number,[],low_bou,up_bou,cheapcon_function);
% x_updata_list=(up_bou'-low_bou').*lhsdesign(sample_number_initial,variable_number)+low_bou';

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

% NFE setting
if isempty(NFE_max)
    if expensive_nonlcon_flag
        NFE_max=50*variable_number;
    else
        NFE_max=20*variable_number;
    end
end

% import data from data library
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);

while ~done
    % step 3
    % updata data library by x_list
    [fval_updata_list,con_updata_list,coneq_updata_list]=dataLibraryUpdata...
        (data_library_name,model_function,x_updata_list);NFE=NFE+size(x_updata_list,1);
    x_list=[x_list;x_updata_list];
    fval_list=[fval_list;fval_updata_list];
    con_list=[con_list;con_updata_list];
    coneq_list=[coneq_list;coneq_updata_list];
    
    % nomalization con average
    fval_max=mean(abs(fval_list),1);
    fval_nomlz_list=fval_list./fval_max*1e3;
    if ~isempty(con_list)
        con_max_list=mean(abs(con_list),1);
        con_nomlz_list=con_list./con_max_list*1e3;
    else
        con_max_list=[];
        con_nomlz_list=[];
    end
    if ~isempty(coneq_list)
        coneq_max_list=mean(abs(coneq_list),1);
        coneq_nomlz_list=coneq_list./coneq_max_list*1e3;
    else
        coneq_max_list=[];
        coneq_nomlz_list=[];
    end
    
    % step 4
    % generate ERBF_QP model use normalization fval
    %     [ERBF_model_fval,ERBF_model_con,ERBF_model_coneq,output_ERBF]=getEnsemleRadialBasisModel...
    %         (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list);
    %     object_function_surrogate=output_ERBF.object_function_surrogate;
    %     nonlcon_function_surrogate=output_ERBF.nonlcon_function_surrogate;
    
    [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output_radialbasis]=getRadialBasisModel...
        (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list);
    object_function_surrogate=output_radialbasis.object_function_surrogate;
    nonlcon_function_surrogate=output_radialbasis.nonlcon_function_surrogate;

    % step 5
    % MSP guideline to obtain x_adapt
    if expensive_nonlcon_flag
        nonlcon_torlance_surrogate=min(max(con_nomlz_list,[],2))*(1-NFE/NFE_max)^2;
    else
        nonlcon_torlance_surrogate=0;
    end
    [x_potential,~,exitflag,~]=findMinMSP...
        (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function_surrogate,...
        cheapcon_function,nonlcon_torlance_surrogate);
    
    if exitflag == -2
        % for expansive nonlcon problem
        object_nonlcon_function_surrogate=@(x) objectNonlconFunctionSurrogate(x,nonlcon_function_surrogate);
        [x_potential,~,exitflag,~]=findMinMSP...
            (object_nonlcon_function_surrogate,variable_number,low_bou,up_bou,[],...
            cheapcon_function,nonlcon_torlance_surrogate);
    end
    
    % check x_potential if exist in data library
    % if not, updata data libraray
    [x_potential,fval_potential,con_potential,coneq_potential,NFE_p]=dataLibraryUpdataProtect...
        (data_library_name,model_function,x_potential',...
        x_list,fval_list,con_list,coneq_list,...
        low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
    x_list=[x_list;x_potential];x_potential=x_potential';
    fval_list=[fval_list;fval_potential];
    if ~isempty(con_list)
        con_list=[con_list;con_potential];con_potential=con_potential';
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_potential];coneq_potential=coneq_potential';
    end
    
    % normalization fval updata
    fval_potential_nomlz=fval_potential/fval_max*1e3;
    fval_nomlz_list=[fval_nomlz_list;fval_potential_nomlz];
    if ~isempty(con_potential)
        con_potential_nomlz=(con_potential'./con_max_list)*1e3;
        con_nomlz_list=[con_nomlz_list;con_potential_nomlz];
    end
    if ~isempty(coneq_potential)
        coneq_potential_nomlz=(coneq_potential'./coneq_max_list)*1e3;
        coneq_nomlz_list=[coneq_nomlz_list;coneq_potential_nomlz];
    end
    
    if DRAW_FIGURE_FLAG && variable_number < 3
        %         interpolationVisualize(ERBF_model_fval,low_bou,up_bou);
        %         line(x_potential(1),x_potential(2),fval_potential_nomlz,'Marker','o','color','r','LineStyle','none')
    end
    
    % step 6
    % find best result to record
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);
    
    if INFORMATION_FLAG
        fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
        fprintf('current x:          %s\n',num2str(x_potential'));
        fprintf('current value:      %f\n',fval_potential);
        fprintf('current violation:  %s  %s\n',num2str(con_potential'),num2str(coneq_potential'));
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
                abs((fval_potential-fval_potential_old)/fval_potential_old) < torlance)
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
    
    % step 7
    if expensive_nonlcon_flag
        % because data prefer getting better
        filter_index_list=[];% filter point list
        feasible_index_list=[];% feasible point list
        
        % generate filter
        for x_index=1:size(x_list,1)
            con_x_nomlz=[];
            if ~isempty(con_nomlz_list)
                con_x_nomlz=max(con_nomlz_list(x_index,:));
            end
            coneq_x_nomlz=[];
            if ~isempty(coneq_nomlz_list)
                coneq_x_nomlz=max(abs(coneq_nomlz_list(x_index,:)));
            end
            total_con_x_nomlz=max([con_x_nomlz;coneq_x_nomlz]);
            
            % only with constraint will add into filter
            if (total_con_x_nomlz > filter_torlance)
                add_filter_flag=1;
                
                filter_index_list_unit=1;
                while filter_index_list_unit <= length(filter_index_list)
                    x_filter_index=filter_index_list(filter_index_list_unit,:);
                    
                    % get h(x) of x and x_filter
                    con_filter_nomlz=[];
                    if ~isempty(con_nomlz_list)
                        con_filter_nomlz=max(con_nomlz_list(x_filter_index,:));
                    end
                    coneq_filter_nomlz=[];
                    if ~isempty(coneq_nomlz_list)
                        coneq_filter_nomlz=max(coneq_nomlz_list(x_filter_index,:));
                    end
                    total_con_filter_nomlz=max([con_filter_nomlz;coneq_filter_nomlz]);
                    
                    % if cannot improve filter, reject it
                    if (fval_nomlz_list(x_index) > fval_nomlz_list(x_filter_index)) && ...
                            (total_con_x_nomlz > total_con_filter_nomlz)
                        add_filter_flag=0;
                        break;
                    end
                    
                    % if better than filter, reject filter
                    if (fval_nomlz_list(x_index) < fval_nomlz_list(x_filter_index)) && ...
                            (total_con_x_nomlz < total_con_filter_nomlz)
                        filter_index_list(filter_index_list_unit)=[];
                        filter_index_list_unit=filter_index_list_unit-1;
                    end
                    
                    filter_index_list_unit=filter_index_list_unit+1;
                end
                % add into filter list if possible
                if add_filter_flag
                    filter_index_list=[filter_index_list;x_index];
                end
            else
                feasible_index_list=[feasible_index_list;x_index];
            end
        end
        %         last_end_index=size(x_list,1);
        if length(feasible_index_list) > 0.2*size(x_list,1)
            filter_torlance=filter_torlance/2;
        end
        
        fval_label=zeros(size(x_list,1),1);
        fval_label(filter_index_list)=1;
        
        % feasible point set label 1
        fval_label(feasible_index_list)=1;
        
        % use filter and train SVM
        kernel_function_SVM=@(x1,x2) exp(-((x1-x2)'*(x1-x2))*scale_SVM);
        [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
            (x_list,fval_label,penalty_SVM,low_bou,up_bou,kernel_function_SVM);
        if DRAW_FIGURE_FLAG && variable_number < 3
            classifySupportVectorMachineVisualization...
                (SVM_model,[13;0],[15;2]);
        end
        
        % get data to obtain clustering center
        x_sup_list=[];
        for x_index=1:sample_number_data
            if  SVM_predict_function(x_data_list(x_index,:)')==1
                x_sup_list=[x_sup_list;x_data_list(x_index,:)];
            end
        end
        
        if isempty(x_sup_list)
            % updata SVM parameter
            if scale_SVM < 1e3
                scale_SVM=scale_SVM*sqrt(10);
            elseif penalty_SVM < 1e3
                penalty_SVM=penalty_SVM*sqrt(10);
            end
            
            % no center found use filter point
            if isempty(feasible_index_list)
                if ~isempty(con_nomlz_list)
                    con_filter_nomlz_list=con_nomlz_list(filter_index_list);
                else
                    con_filter_nomlz_list=[];
                end
                if ~isempty(coneq_nomlz_list)
                    coneq_filter_nomlz_list=coneq_nomlz_list(filter_index_list);
                else
                    coneq_filter_nomlz_list=[];
                end
                max_totalcon_list=max([con_filter_nomlz_list,coneq_filter_nomlz_list],[],2);
                [~,filter_min_index]=min(max_totalcon_list);
                x_center=x_list(filter_index_list(filter_min_index),:)';
            else
                [~,min_fval_index]=min(fval_list(feasible_index_list));
                x_center=x_list(feasible_index_list(min_fval_index),:)';
            end
        end
    else
        % interset sampling
        delta_add=0.1;
        fval_thresh=min(fval_list)+(eta-delta_add)*(max(fval_list)-min(fval_list));
        x_sup_list=[];
        
        while size(x_sup_list,1) < 5*variable_number
            % step 7-1
            % classify exist data
            fval_thresh=fval_thresh+delta_add*(max(fval_list)-min(fval_list));
            fval_label=zeros(size(x_list,1),1);
            for x_index=1:size(x_list,1)
                if fval_list(x_index) <= fval_thresh
                    fval_label(x_index)=1;
                else
                    fval_label(x_index)=0;
                end
            end
            
            if sum(fval_label) == size(x_list,1)
                delta_add=delta_add/2;
                fval_thresh=min(fval_list)+delta_add*(max(fval_list)-min(fval_list));
            else
                kernel_function_SVM=@(x1,x2) exp(-((x1-x2)'*(x1-x2))*scale_SVM);
                % step 7-2
                % get a large number of x point, use SVM to predict x point
                [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
                    (x_list,fval_label,penalty_SVM,low_bou,up_bou,kernel_function_SVM);
                if DRAW_FIGURE_FLAG && variable_number < 3
                    classifySupportVectorMachineVisualization...
                        (SVM_model,low_bou,up_bou);
                end
                for x_index=1:sample_number_data
                    if  SVM_predict_function(x_data_list(x_index,:)')==1
                        x_sup_list=[x_sup_list;x_data_list(x_index,:)];
                    end
                end
            end
        end
    end
    
    % step 7-3
    % calculate clustering center
    if ~isempty(x_sup_list)
        FC_model=classifyFuzzyClusteringFeatureSpace...
            (x_sup_list,1,low_bou,up_bou,m);
        x_center=FC_model.center_list';
    end
    
    % updata ISR
    x_potential_nomlz=(x_potential-low_bou)./(up_bou-low_bou);
    x_center_nomlz=(x_center-low_bou)./(up_bou-low_bou);
    bou_range_nomlz=eta*norm(x_potential_nomlz-x_center_nomlz,2);
    if bou_range_nomlz < 1e-2
        bou_range_nomlz=1e-2;
    end
    bou_range=bou_range_nomlz.*(up_bou-low_bou);
    low_bou_ISR=x_potential-bou_range;
    low_bou_ISR=max(low_bou_ISR,low_bou);
    up_bou_ISR=x_potential+bou_range;
    up_bou_ISR=min(up_bou_ISR,up_bou);
    
    if DRAW_FIGURE_FLAG && variable_number < 3
        bou_line=[low_bou_ISR,[low_bou_ISR(1);up_bou_ISR(2)],up_bou_ISR,[up_bou_ISR(1);low_bou_ISR(2)],low_bou_ISR];
        line(bou_line(1,:),bou_line(2,:));
        line(x_potential(1),x_potential(2),'Marker','x')
    end
    
    % sampling in ISR
    [x_list_exist,~,~,~]=dataLibraryLoad...
        (data_library_name,low_bou_ISR,up_bou_ISR);
    [~,x_updata_list,~]=getLatinHypercube...
        (sample_number_iteration+size(x_list_exist,1),variable_number,x_list_exist,...
        low_bou_ISR,up_bou_ISR,cheapcon_function);
    %     x_updata_list=(up_bou_ISR'-low_bou_ISR').*lhsdesign(sample_number_iteration,variable_number)+low_bou_ISR';
    
    x_potential_old=x_potential;
    fval_potential_old=fval_potential;
    fval_best_old=fval_best;
end
result_x_best=result_x_best(1:iteration-1,:);
result_fval_best=result_fval_best(1:iteration-1);

output.result_x_best=result_x_best;
output.result_fval_best=result_fval_best;

    function [fval,con,coneq]=modelFunction...
            (x,object_function,nonlcon_function)
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
    function fval=objectNonlconFunctionSurrogate(x,nonlcon_function_surrogate)
        [con__,coneq__]=nonlcon_function_surrogate(x);
        fval=0;
        if ~isempty(con__)
            fval=fval+sum(max(con__,0).^2);
        end
        if ~isempty(coneq__)
            fval=fval+sum(max(con__,0).^2);
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

%% FCM
function FC_model=classifyFuzzyClusteringFeatureSpace...
    (X,classify_number,low_bou,up_bou,m,kernal_function)
% get fuzzy cluster model with feature space
% kernal function recommend kernal_function=@(sq) exp(-sq/2*1000);
% X is x_number x variable_number matrix
% center_list is classify_number x variable_number matrix
%
if nargin < 6
    kernal_function=[];
    if nargin < 4
        up_bou=[];
        if nargin < 3
            low_bou=[];
        end
    end
end
iteration_max=100;
torlance=1e-3;

[x_number,variable_number]=size(X);

if isempty(kernal_function)
    kernal_function=@(sq) exp(-sq/2*1000);
end

% nomalization data
if isempty(low_bou)
    low_bou=min(X,[],1)';
end
if isempty(up_bou)
    up_bou=max(X,[],1)';
end
X_nomlz=(X-low_bou')./(up_bou'-low_bou');

% if x_number equal 1, clustering cannot done
if x_number==1
    FC_model.X=X;
    FC_model.X_normalize=X_nomlz;
    FC_model.center_list=X;
    FC_model.fval_loss_list=[];
    return;
end

U=zeros(classify_number,x_number);
center_list=rand(classify_number,variable_number)*0.5;
iteration=0;
done=0;
fval_loss_list=zeros(iteration_max,1);

% get X_center_dis_sq
X_center_dis_sq=zeros(classify_number,x_number);
for classify_index=1:classify_number
    for x_index=1:x_number
        X_center_dis_sq(classify_index,x_index)=...
            getSq((X_nomlz(x_index,:)-center_list(classify_index,:)));
    end
end

while ~done
    % updata classify matrix U
    for classify_index=1:classify_number
        for x_index=1:x_number
            U(classify_index,x_index)=...
                1/sum(((2-2*kernal_function(X_center_dis_sq(classify_index,x_index)))./...
                (2-2*kernal_function(X_center_dis_sq(:,x_index)))).^(1/(m-1)));
        end
    end
    
    % updata center_list
    center_list_old=center_list;
    for classify_index=1:classify_number
        center_list(classify_index,:)=...
            sum((U(classify_index,:)').^m.*X_nomlz,1)./...
            sum((U(classify_index,:)').^m,1);
    end
    
    % updata X_center_dis_sq
    X_center_dis_sq=zeros(classify_number,x_number);
    for classify_index=1:classify_number
        for x_index=1:x_number
            X_center_dis_sq(classify_index,x_index)=...
                getSq((X_nomlz(x_index,:)-center_list(classify_index,:)));
        end
    end
    
    % forced interrupt
    if iteration > iteration_max
        done=1;
    end
    
    % convergence judgment
    if sum(sum(center_list_old-center_list).^2)<torlance
        done=1;
    end
    
    iteration=iteration+1;
    fval_loss_list(iteration)=sum(sum(U.^m.*(2-2*kernal_function(X_center_dis_sq))));
end
fval_loss_list(iteration+1:end)=[];
center_list=center_list.*(up_bou'-low_bou')+low_bou';

FC_model.X=X;
FC_model.X_normalize=X_nomlz;
FC_model.center_list=center_list;
FC_model.fval_loss_list=fval_loss_list;

    function sq=getSq(dx)
        % dx is 1 x variable_number matrix
        %
        sq=dx*dx';
    end
end

%% multi-fidelity gaussian process classifier
function [predict_function,CGPMF_model]=classifyGaussProcessMultiFidelity...
    (XHF,ClassHF,XLF,ClassLF,low_bou,up_bou)
% generate gauss classifier model
% version 6,this version is assembly of gpml-3.6 EP method
% X, XHF, XLF is x_number x variable_number matirx
% Class, ClassH, ClassLF is x_number x 1 matrix
% low_bou, up_bou is 1 x variable_number matrix
%
% abbreviation:
% pred: predicted,nomlz: normalization,num: number
% var: variance
%
X={XHF,XLF};
Class={ClassHF,ClassLF};
if nargin < 4 || isempty(up_bou)
    up_bou=max([XHF;XLF],[],1);
end
if nargin < 3 || isempty(low_bou)
    low_bou=min([XHF;XLF],[],1);
end

[~,variable_number]=size(X);

% normalization data
X_nomlz={(XHF-low_bou)./(up_bou-low_bou),...
    (XLF-low_bou)./(up_bou-low_bou)};

object_function=@(x) objectFunctionGPC(x,{@infEP},{@meanConst},{@calCovMF},{@likErf},X_nomlz,Class);
x=zeros(1,2*variable_number+2);
low_bou_hyp=log([ones(1,2*variable_number+2)*0.1]);
up_bou_hyp=log([ones(1,2*variable_number+2)*10]);
% up_bou_hyp(end)=0;
x=fmincon(object_function,x,[],[],[],[],low_bou_hyp,up_bou_hyp,[],...
    optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',10));

hyp.mean=x(1);
hyp.cov=x(2:end);
hyp.lik=[];
post=infEP(hyp,{@meanConst},{@calCovMF},{@likErf},X_nomlz,Class);
predict_function=@(x_pred) classifyGaussPredictor...
    (x_pred,hyp,{@meanConst},{@calCovMF},{@likErf},post,X_nomlz,low_bou,up_bou);

% output model
CGPMF_model.X=X;
CGPMF_model.Class=Class;
CGPMF_model.X_nomlz=X_nomlz;
CGPMF_model.low_bou=low_bou;
CGPMF_model.up_bou=up_bou;
CGPMF_model.predict_function=predict_function;
CGPMF_model.hyp=hyp;
CGPMF_model.post=post;

    function [fval,gradient]=objectFunctionGPC(x,inf,mean,cov,lik,X,Y)
        hyp_iter.mean=x(1);
        hyp_iter.cov=x(2:end);
        hyp_iter.lik=[];

        if nargout < 2
            [~,nlZ] = feval(inf{:},hyp_iter,mean,cov,lik,X,Y);
            fval=nlZ;
        elseif nargout < 3
            [~,nlZ,dnlZ]=feval(inf{:},hyp_iter,mean,cov,lik,X,Y);
            fval=nlZ;
            gradient=[dnlZ.mean,dnlZ.cov];
        end
    end
    function [class,possibility,miu_pre,var_pre]=classifyGaussPredictor...
            (x_pred,hyp,mean,cov,lik,post,X,low_bou,up_bou)
        % predict function
        %
        x_pred_nomlz=(x_pred-low_bou)./(up_bou-low_bou);
        pred_num=size(x_pred_nomlz,1);
        ys=ones(pred_num,1);

        alpha = post.alpha; L = post.L; sW = post.sW;
        %verify whether L contains valid Cholesky decomposition or something different
        Lchol = isnumeric(L) && all(all(tril(L,-1)==0)&diag(L)'>0&isreal(diag(L))');
        ns = size(x_pred_nomlz,1);                                       % number of data points
        nperbatch = 1000;                       % number of data points per mini batch
        nact = 0;                       % number of already processed test data points
        ymu = zeros(ns,1); ys2 = ymu; miu_pre = ymu; var_pre = ymu; possibility = ymu;   % allocate mem
        while nact<ns               % process minibatches of test cases to save memory
            id = (nact+1):min(nact+nperbatch,ns);               % data points to process
            kss = feval(cov{:},hyp.cov,x_pred_nomlz(id,:),'diag');              % self-variance
            Ks = feval(cov{:},hyp.cov,X,x_pred_nomlz(id,:));        % avoid computation
            ms = feval(mean{:},hyp.mean,x_pred_nomlz(id,:));
            N = size(alpha,2);  % number of alphas (usually 1; more in case of sampling)
            Fmu = repmat(ms,1,N) + Ks'*full(alpha);        % conditional mean fs|f
            miu_pre(id) = sum(Fmu,2)/N;                                   % predictive means
            if Lchol    % L contains chol decomp => use Cholesky parameters (alpha,sW,L)
                V  = L'\(repmat(sW,1,length(id)).*Ks);
                var_pre(id) = kss - sum(V.*V,1)';                       % predictive variances
            else                % L is not triangular => use alternative parametrisation
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

        possibility=exp(possibility);
        class=ones(pred_num,1);
        index_list=find(possibility < 0.5);
        class(index_list)=-1;
    end
end

function [K,dK_dvar]=calCovMF(cov,X,Z)
% obtain covariance of x
%
if iscell(X)
    variable_number=size(X{1},2);
else
    variable_number=size(X,2);
end

lenH=exp(cov(1:variable_number));
lenL=exp(cov(variable_number+(1:variable_number)));
rho=exp(cov(end));

if nargin > 2 && nargout < 2 && ~isempty(Z)
    if strcmp(Z,'diag')
        K=rho*rho*1+1;
        return
    end
end

XHF=X{1};
XLF=X{2};
[xH_num,variable_number]=size(XHF);
[xL_num,~]=size(XLF);
x_num=xH_num+xL_num;
X=[XHF;XLF];

% predict
if nargin > 2 && nargout < 2 && ~isempty(Z)
    [z_num,variable_number]=size(Z);
    % initializate square of X inner distance
    sq_dis=zeros(x_num,z_num,variable_number);
    for len_index=1:variable_number
        sq_dis(:,:,len_index)=(X(:,len_index)-Z(:,len_index)').^2;
    end

    % exp of x__x with H
    exp_disH=zeros(xH_num,z_num);
    for len_index=1:variable_number
        exp_disH=exp_disH+...
            sq_dis(1:xH_num,:,len_index)/2/lenH(len_index)^2;
    end
    exp_disH=exp(-exp_disH);

    % exp of x__x with L
    exp_disL=zeros(x_num,z_num);
    for len_index=1:variable_number
        exp_disL=exp_disL+...
            sq_dis(1:x_num,:,len_index)/2/lenL(len_index)^2;
    end
    exp_disL=exp(-exp_disL);

    K=exp_disL;
    K(1:xH_num,:)=rho*rho*K(1:xH_num,:)+exp_disH;
    K(xH_num+1:end,:)=rho*K(xH_num+1:end,:);
else
    % initializate square of X inner distance
    sq_dis=zeros(x_num,x_num,variable_number);
    for len_index=1:variable_number
        sq_dis(:,:,len_index)=(X(:,len_index)-X(:,len_index)').^2;
    end

    % exp of x__x with H
    exp_disH=zeros(xH_num);
    for len_index=1:variable_number
        exp_disH=exp_disH+...
            sq_dis(1:xH_num,1:xH_num,len_index)/2/lenH(len_index)^2;
    end
    exp_disH=exp(-exp_disH);
    KH=exp_disH;

    % exp of x__x with L
    exp_disL=zeros(x_num);
    for len_index=1:variable_number
        exp_disL=exp_disL+...
            sq_dis(1:end,1:end,len_index)/2/lenL(len_index)^2;
    end
    exp_disL=exp(-exp_disL);
    % times rho: HH to rho2, HL to rho, LL to 1
    rho_exp_disL=exp_disL;
    rho_exp_disL(1:xH_num,1:xH_num)=...
        (rho*rho)*exp_disL(1:xH_num,1:xH_num);
    rho_exp_disL(1:xH_num,(xH_num+1):end)=...
        rho*exp_disL(1:xH_num,(xH_num+1):end);
    rho_exp_disL((xH_num+1):end,1:xH_num)=...
        rho_exp_disL(1:xH_num,(xH_num+1):end)';

    KL=rho_exp_disL;
    K=KL;
    K(1:xH_num,1:xH_num)=K(1:xH_num,1:xH_num)+KH;

    if nargout >= 2
        dK_dvar=cell(1,2*variable_number+1);

        % len H
        for len_index=1:variable_number
            dK_dlenH=zeros(x_num);
            dK_dlenH(1:xH_num,1:xH_num)=KH.*...
                sq_dis(1:xH_num,1:xH_num,len_index)/lenH(len_index)^2;
            dK_dvar{len_index}=dK_dlenH;
        end

        % len L
        for len_index=1:variable_number
            dK_dlenL=KL.*sq_dis(:,:,len_index)/lenL(len_index)^2;
            dK_dvar{(variable_number)+len_index}=dK_dlenL;
        end

        % rho
        dK_drho=zeros(x_num);
        dK_drho(1:xH_num,1:xH_num)=...
            2*rho*rho*exp_disL(1:xH_num,1:xH_num);
        dK_drho((xH_num+1):end,1:end)=...
            rho*exp_disL((xH_num+1):end,1:end);
        dK_drho(1:end,(xH_num+1):end)=...
            dK_drho((xH_num+1):end,1:end)';
        dK_dvar{end}=dK_drho;
    end
end

end

function [post,nlZ,dnlZ] = infEP(hyp,mean,cov,lik,xF,yF)
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
if isnumeric(cov), K = cov;                    % use provided covariance matrix
else K = feval(cov{:}, hyp.cov, xF); end       % evaluate the covariance matrix

x=[xF{1};xF{2}];
y=[yF{1};yF{2}];
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

if ~isreal(tau_ni)
disp('?');
end

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

if sweep == max_sweep && abs(nlZ-nlZ_old) > tol
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
    [K,dK] = feval(cov{:},hyp.cov,xF,[]);
    for i=1:length(hyp.cov)
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
L = chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
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

if ~isreal(Sigma(1))
disp('?');
end
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
if nargin==2
    A = c*ones(size(x,1),1);                                       % evaluate mean
else
    if i==1
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
if nargin>1,y = sign(y); y(y==0) = 1; else y = 1; end % allow only +/- 1 values
if numel(y)==0,y = 1; end

if nargin<5                              % prediction mode if inf is not present
    y = y.*ones(size(mu));                                       % make y a vector
    s2zero = 1; if nargin>3&&numel(s2)>0&&norm(s2)>eps,s2zero = 0; end  % s2==0 ?
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
                if numel(y)>0,z = z.*y; end
                if nargout<=1,lZ = logphi(z);                         % log part function
                else          [lZ,n_p] = logphi(z); end
                if nargout>1
                    if numel(y)==0,y=1; end
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
f = 0; for i=1:14,f = lp0.*(c(i)+f); end,lp(id1) = -2*f-log(2);
id2 = z<-11.3137;                                    % second case: very small
r = [ 1.2753666447299659525; 5.019049726784267463450;
    6.1602098531096305441; 7.409740605964741794425;
    2.9788656263939928886 ];
q = [ 2.260528520767326969592;  9.3960340162350541504;
    12.048951927855129036034; 17.081440747466004316;
    9.608965327192787870698;  3.3690752069827527677 ];
num = 0.5641895835477550741; for i=1:5,num = -z(id2).*num/sqrt(2) + r(i); end
den = 1.0;                   for i=1:6,den = -z(id2).*den/sqrt(2) + q(i); end
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

%% multi-fidelity surrogate model
function [MK_model_fval,MK_model_con,MK_model_coneq,output]=getHieraKrigingModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create radialbasis model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
MK_model_fval=interpolationRadialBasisPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    MK_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'radialbasis_matrix',[],'inv_radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    MK_model_con=repmat(MK_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        MK_model_con(con_index)=interpolationRadialBasisPreModel...
            (x_list,con_list(:,con_index));
    end
else
    MK_model_con=[];
end

if ~isempty(coneq_list)
    MK_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
        'radialbasis_matrix',[],'inv_radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    MK_model_coneq=repmat(MK_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        MK_model_coneq(coneq_index)=interpolationRadialBasisPreModel...
            (x_list,coneq_list(:,coneq_index));
    end
else
    MK_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,MK_model_fval);
if isempty(MK_model_con) && isempty(MK_model_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,MK_model_con,MK_model_coneq);
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

function [predict_function,HK_model]=interpHieraKrigingPreModel...
    (XHF,YHF,XLF,YLF,hyp)
% construct Hierarchical Kriging version 1
% XHF, YHF are x_HF_number x variable_number matrix
% XLF, YLF are x_LF_number x variable_number matrix
% aver_X,stdD_X is 1 x x_HF_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
% hyp: hyp_HF, hyp_LF
% notice theta=exp(hyp)
%
% input:
% XHF, YHF, XLF, YLF, hyp(hyp_HF, hyp_LF)
%
% output:
% kriging model (predict_function,
% XHF, YHF, XLF, YLF, base_function_list)
%
% reference: [1] HAN Z-H, GöRTZ S. Hierarchical Kriging Model for
% Variable-Fidelity Surrogate Modeling [J]. AIAA Journal, 2012, 50(9):
% 1885-96.
%
% Copyright 2023.2 Adel
%
X=[XHF;XLF];
Y=[YHF;YLF];
[x_number,variable_number]=size(X);
x_HF_number=size(XHF,1);
x_LF_number=size(XLF,1);
if nargin < 5
    hyp=ones(1,2*variable_number);
end
hyp_LF=hyp((variable_number+1):end);
hyp_HF=hyp(1:variable_number);

% normalize data
aver_X=mean(X);
stdD_X=std(X);
aver_Y=mean(Y);
stdD_Y=std(Y);
index__=find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__)=1; end
index__=find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__)=1; end
XHF_nomlz=(XHF-aver_X)./stdD_X;
YHF_nomlz=(YHF-aver_Y)./stdD_Y;
XLF_nomlz=(XLF-aver_X)./stdD_X;
YLF_nomlz=(YLF-aver_Y)./stdD_Y;

% first step
% construct low fidelity model

% initial X_dis_sq
XHF_dis_sq=zeros(x_HF_number,x_HF_number,variable_number);
for variable_index=1:variable_number
    XHF_dis_sq(:,:,variable_index)=...
        (XHF_nomlz(:,variable_index)-XHF_nomlz(:,variable_index)').^2;
end
XLF_dis_sq=zeros(x_LF_number,x_LF_number,variable_number);
for variable_index=1:variable_number
    XLF_dis_sq(:,:,variable_index)=...
        (XLF_nomlz(:,variable_index)-XLF_nomlz(:,variable_index)').^2;
end

% regression function define
% notice reg_function process no normalization data
% reg_function=@(X) regZero(X);
reg_function=@(X) regLinear(X);

% calculate reg
fval_reg_nomlz_LF=(reg_function(XLF)-aver_Y)./stdD_Y;

% optimal to get hyperparameter
low_bou_hyp=-4*ones(1,variable_number);
up_bou_hyp=4*ones(1,variable_number);
fmincon_option=optimoptions('fmincon','Display','none',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',true);
object_function_hyp=@(hyp) objectNLLKriging...
    (XLF_dis_sq,YLF_nomlz,x_LF_number,variable_number,hyp,fval_reg_nomlz_LF);

% [fval,gradient]=object_function_hyp(hyp_LF)
% [~,gradient_differ]=differ(object_function_hyp,hyp_LF)

hyp_LF=fmincon...
    (object_function_hyp,hyp_LF,[],[],[],[],low_bou_hyp,up_bou_hyp,[],fmincon_option);

% get parameter
[cov_LF,inv_cov_LF,beta_LF,sigma_sq_LF]=interpKriging...
    (XLF_dis_sq,YLF_nomlz,x_LF_number,variable_number,exp(hyp_LF),fval_reg_nomlz_LF);
gama_LF=inv_cov_LF*(YLF_nomlz-fval_reg_nomlz_LF*beta_LF);
FTRF_LF=fval_reg_nomlz_LF'*inv_cov_LF*fval_reg_nomlz_LF;

% initialization predict function
predict_function_LF=@(X_predict) interpKrigingPredictor...
    (X_predict,XLF_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_LF_number,variable_number,exp(hyp_LF),beta_LF,gama_LF,sigma_sq_LF,...
    inv_cov_LF,fval_reg_nomlz_LF,FTRF_LF,reg_function);

% second step
% construct hierarchical model

% evaluate low fidelty predict value in high fidelity point as base fval
reg_function=@(X) predict_function_LF(X);

% calculate reg
fval_reg_nomlz_HF=(reg_function(XHF)-aver_Y)./stdD_Y;

% optimal to get hyperparameter
object_function_hyp=@(hyp) objectNLLKriging...
    (XHF_dis_sq,YHF_nomlz,x_HF_number,variable_number,hyp,fval_reg_nomlz_HF);

% [fval,gradient]=object_function_hyp(hyp_HF)
% [~,gradient_differ]=differ(object_function_hyp,hyp_HF)

% drawFunction(object_function_hyp,low_bou_hyp,up_bou_hyp);

hyp_HF=fmincon...
    (object_function_hyp,hyp_HF,[],[],[],[],low_bou_hyp,up_bou_hyp,[],fmincon_option);

% calculate covariance and other parameter
[cov_HF,inv_cov_HF,beta_HF,sigma_sq_HF]=interpKriging...
    (XHF_dis_sq,YHF_nomlz,x_HF_number,variable_number,exp(hyp_HF),fval_reg_nomlz_HF);
gama_HF=inv_cov_HF*(YHF_nomlz-fval_reg_nomlz_HF*beta_HF);
FTRF_HF=fval_reg_nomlz_HF'*inv_cov_HF*fval_reg_nomlz_HF;

% initialization predict function
predict_function=@(X_predict) interpKrigingPredictor...
    (X_predict,XHF_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_HF_number,variable_number,exp(hyp_HF),beta_HF,gama_HF,sigma_sq_HF,...
    inv_cov_HF,fval_reg_nomlz_HF,FTRF_HF,reg_function);

HK_model.X=XHF;
HK_model.Y=YHF;
HK_model.XHF=XHF;
HK_model.YHF=YHF;
HK_model.XLF=XLF;
HK_model.YLF=YLF;
HK_model.cov_LF=cov_LF;
HK_model.inv_cov_LF=inv_cov_LF;
HK_model.cov_HF=cov_HF;
HK_model.inv_cov_HF=inv_cov_HF;

hyp=[hyp_HF,hyp_LF];
HK_model.hyp=hyp;
HK_model.aver_X=aver_X;
HK_model.stdD_X=stdD_X;
HK_model.aver_Y=aver_Y;
HK_model.stdD_Y=stdD_Y;

HK_model.predict_function=predict_function;

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
        cov=exp(-cov/vari_num)+eye(x_num)*1e-6;

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
            inv_cov,fval_reg,FTRF,reg_function)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num,~]=size(X_pred);
        fval_reg_pred=reg_function(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;
        fval_reg_pred_nomlz=(fval_reg_pred-aver_Y)./stdD_Y;
        
        % predict covariance
        predict_cov=zeros(x_num,x_pred_num);
        for vari_index=1:vari_num
            predict_cov=predict_cov+...
                (X_nomlz(:,vari_index)-X_pred_nomlz(:,vari_index)').^2*theta(vari_index);
        end
        predict_cov=exp(-predict_cov);

        % predict base fval
        
        Y_pred=fval_reg_pred_nomlz*beta+predict_cov'*gama;
        
        % predict variance
        u__=fval_reg'*inv_cov*predict_cov-fval_reg_pred_nomlz';
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

%% data library
function [fval_list,con_list,coneq_list]=dataLibraryUpdata...
    (data_library_name,model_function,x_list)
% updata data library
% updata format:
% variable_number, fval_number, con_number, coneq_number
% x, fval, con, coneq
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
file_optimalSurrogate_output=fopen(data_library_name,'a');
file_result=fopen('result_total.txt','a');

% updata format:
% variable_number, fval_number, con_number, coneq_number
% x, fval, con, coneq
for x_index=1:x_number
    x=x_list(x_index,:);
    [fval,con,coneq]=model_function(x);
    fval=fval(:);
    con=con(:);
    coneq=coneq(:);
    fval_list=[fval_list;fval'];
    con_list=[con_list;con'];
    coneq_list=[coneq_list;coneq'];
    
    % write data to txt_optimalSurrogateSADEKTS
    fprintf(file_optimalSurrogate_output,'%d ',variable_number);
    fprintf(file_optimalSurrogate_output,'%d ',length(fval));
    fprintf(file_optimalSurrogate_output,'%d ',length(con));
    fprintf(file_optimalSurrogate_output,'%d ',length(coneq));
    
    fprintf(file_optimalSurrogate_output,x_format,x);
    fval_format=repmat(fval_format_base,1,length(fval));
    fprintf(file_optimalSurrogate_output,fval_format,fval);
    fval_format=repmat(fval_format_base,1,length(con));
    fprintf(file_optimalSurrogate_output,fval_format,con);
    fval_format=repmat(fval_format_base,1,length(coneq));
    fprintf(file_optimalSurrogate_output,fval_format,coneq);
    fprintf(file_optimalSurrogate_output,'\n');
    
    % write data to txt_result
    fprintf(file_result,'%d ',variable_number);
    fprintf(file_result,'%d ',length(fval));
    fprintf(file_result,'%d ',length(con));
    fprintf(file_result,'%d ',length(coneq));
    
    fprintf(file_result,x_format,x);
    fval_format=repmat(fval_format_base,1,length(fval));
    fprintf(file_result,fval_format,fval);
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
% updata format:
% variable_number, fval_number, con_number, coneq_number
% x, fval, con, coneq
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

% updata format:
% variable_number, fval_number, con_number, coneq_number
% x, fval, con, coneq
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
            fval_number=data(2);
            con_number=data(3);
            coneq_number=data(4);
            
            base=5;
            x=data(base:base+variable_number-1);
            judge=sum(x < low_bou)+sum(x > up_bou);
            if ~judge
                x_list=[x_list;x];
                base=base+variable_number;
                fval_list=[fval_list;data(base:base+fval_number-1)];
                base=base+fval_number;
                con=data(base:base+con_number-1);
                if ~isempty(con)
                    con_list=[con_list;con];
                end
                base=base+con_number;
                coneq=data(base:base+coneq_number-1);
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
% election sequential method is used(sample and iteration)
%
% sample number is total point in area
% default low_bou is 0, up_bou is 1, cheapcon_function is []
% low_bou and up_bou is colume vector
% X_exist, X, supply_X is x_number x variable_number matrix
% X_exist should meet bou
%
% reference:
% [1]LONG T, LI X, SHI R, et al., Gradient-Free Trust-Region-Based
% Adaptive Response Surface Method for Expensive Aircraft Optimization[J].
% AIAA Journal, 2018, 56(2): 862-73.
%
% Copyright 2022 Adel
%
if nargin < 6
    cheapcon_function=[];
    if nargin < 5
        if nargin < 3
            if nargin < 2
                error('getLatinHypercube: lack variable_number');
            end
        end
        low_bou=zeros(variable_number,1);
        up_bou=ones(variable_number,1);
    end
end

% check x_exist_list if meet boundary
if ~isempty(X_exist)
    if size(X_exist,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
    index=find(X_exist < low_bou);
    index=[index,find(X_exist > up_bou)];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list range error');
    end
    X_exist_nomlz=(X_exist-low_bou)./(up_bou-low_bou);
else
    X_exist_nomlz=[];
end

if sample_number <= 0
    X=[];
    X_new=[];
    distance_min_nomlz=[];
    return;
end

iteration_max=1000*variable_number;
x_new_number=sample_number-size(X_exist,1);
if x_new_number <= 0
    X=X_exist;
    X_new=[];
    distance_min_nomlz=getMinDistance(X_exist_nomlz);
    return;
end

% get quasi-feasible point
x_initial_number=100*x_new_number;
if ~isempty(cheapcon_function)
    X_supply_quasi_nomlz=[];
    % check if have enough X_supply_nomlz
    while size(X_supply_quasi_nomlz,1) < 100*x_new_number
        X_supply_initial_nomlz=rand(x_initial_number,variable_number);
        x_index=1;
        while x_index <= size(X_supply_initial_nomlz,1)
            x_supply=X_supply_initial_nomlz(x_index,:).*(up_bou-low_bou)+low_bou;
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
distance_min_nomlz=0;
X_new_nomlz=[];
while iteration <= iteration_max
    % random select x_new_number X to X_trial_nomlz
    x_select_index=randperm(x_supply_quasi_number,x_new_number);
    
    % get distance min itertion X_
    distance_min_iteration=getMinDistanceIter...
        (X_supply_quasi_nomlz(x_select_index,:),X_exist_nomlz);
    
    % if distance_min_iteration is large than last time
    if distance_min_iteration > distance_min_nomlz
        distance_min_nomlz=distance_min_iteration;
        X_new_nomlz=X_supply_quasi_nomlz(x_select_index,:);
    end
    
    iteration=iteration+1;
end
X_new=X_new_nomlz.*(up_bou-low_bou)+low_bou;
X=[X_new;X_exist];
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);

    function distance_min__=getMinDistance(x_list__)
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
