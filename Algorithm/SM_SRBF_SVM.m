clc;
clear;
close all hidden;

benchmark=BenchmarkFunction();

% variable_number=2;
% object_function=@(x) benchmark.singleGPObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-2,-2];
% up_bou=[2,2];
% nonlcon_function=[];
% cheapcon_function=[];

% variable_number=2;
% object_function=@(x) benchmark.singlePKObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-3,-3];
% up_bou=[3,3];
% nonlcon_function=[];
% cheapcon_function=[];

% variable_number=6;
% object_function=@(x) benchmark.singleHNObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=zeros(1,variable_number);
% up_bou=ones(1,variable_number);
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

% variable_number=20;
% object_function=@(x) benchmark.singleDP20Object(x);
% object_function_LF=@(x) benchmark.singleDP20ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=ones(1,variable_number)*-30;
% up_bou=ones(1,variable_number)*30;
% nonlcon_function=[];
% nonlcon_function_LF=[];
% cheapcon_function=[];

% variable_number=20;
% object_function=@(x) benchmark.singleEP20Object(x);
% object_function_LF=@(x) benchmark.singleEP20ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=ones(1,variable_number)*-30;
% up_bou=ones(1,variable_number)*30;
% nonlcon_function=[];
% nonlcon_function_LF=[];
% cheapcon_function=[];

% variable_number=2;
% object_function=@(x) benchmark.singleG06Object(x);
% object_function_LF=@(x) benchmark.singleG06ObjectLow(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[13,0];
% up_bou=[100,100];
% nonlcon_function=@(x) benchmark.singleG06Nonlcon(x);
% nonlcon_function_LF=@(x) benchmark.singleG06NonlconLow(x);
% cheapcon_function=[];
% model_function=[];

variable_number=4;
object_function=@(x) benchmark.singlePVD4Object(x);
object_function_low=@(x) benchmark.singlePVD4ObjecttLow(x);
A=[-1,0,0.0193,0;
    0,-1,0.00954,0;];
B=[0;0];
Aeq=[];
Beq=[];
low_bou=[0,0,0,0];
up_bou=[1,1,50,240];
nonlcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq,@(x) benchmark.singlePVD4Nonlcon(x));
cheapcon_function=[];
model_function=[];

% variable_number=13;
% object_function=@benchmark.singleG01Object;
% object_function_low=@benchmark.singleG01ObjectLow;
% A=[ 2   2   0   0   0   0   0   0   0   1   1   0   0;
%     2   0   2   0   0   0   0   0   0   1   0   1   0;
%     0   2   2   0   0   0   0   0   0   0   1   1   0;
%     -8  0   0   0   0   0   0   0   0   1   0   0   0;
%     0   -8  0   0   0   0   0   0   0   0   1   0   0;
%     0   0   0   -2  -1  0   0   0   0   1   0   0   0;
%     0   0   0   0   0   -2  -1  0   0   0   1   0   0;
%     0   0   0   0   0   0   0   -2  -1  0   0   1   0;
%     ];
% B=[10;10;10;0;0;0;0;0];
% Aeq=[];
% Beq=[];
% low_bou=zeros(1,13);
% up_bou=ones(1,13);
% up_bou(10:12)=100;
% nonlcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq,[]);
% nonlcon_function_LF=@(x) cheapconFunction(x,A,B,Aeq,Beq,[]);
% cheapcon_function=[];

% x_initial=rand(1,variable_number).*(up_bou-low_bou)+low_bou;
% fmincon_option=optimoptions('fmincon','Algorithm','sqp');
% [x_best,fval_best,~,output,lambda,grad,hessian]=fmincon(object_function,x_initial,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,fmincon_option)

data_library_name='optimal_data_library.txt';
delete(data_library_name);
delete('result_total.txt');
[x_best,fval_best,NFE,output]=optimalSurrogateSRBFSVM...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,[],50)

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

%% main
function [x_best,fval_best,NFE,output]=optimalSurrogateSRBFSVM...
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
% referance: [1] SHI R, LIU L, LONG T, et al. Sequential Radial Basis
% Function Using Support Vector Machine for Expensive Design Optimization
% [J]. AIAA Journal, 2017, 55(1): 214-27.
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
CONVERGENCE_JUDGMENT_FLAG=1; % whether judgment convergence

if isempty(iteration_max)
    iteration_max=100;
end

% hyper parameter
sample_number_initial=min((variable_number+1)*(variable_number+2)/2,5*variable_number);
sample_number_iteration=variable_number;
sample_number_data=10*sample_number_initial;
eta=1/variable_number; % space decrease coefficient

penalty_SVM=100;
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
    (up_bou-low_bou)+low_bou;

% if do not input model_function, generate model_function
if isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end

iteration=iteration+1;

% step 2
% use latin hypercube method to get initial sample x_list
% [~,x_updata_list,~]=getLatinHypercube...
%     (sample_number_initial,variable_number,[],low_bou,up_bou,cheapcon_function);
x_updata_list=lhsdesign(sample_number_initial,variable_number).*(up_bou-low_bou)+low_bou;

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
    if ~isempty(con_updata_list)
        con_list=[con_list;con_updata_list];
    end
    if ~isempty(coneq_updata_list)
        coneq_list=[coneq_list;coneq_updata_list];
    end

    % nomalization all data
    fval_max=max(abs(fval_list),[],1);
    fval_nomlz_list=fval_list./fval_max*10;
    if ~isempty(con_list)
        con_max_list=max(abs(con_list),[],1);
        con_nomlz_list=con_list./con_max_list*10;
    else
        con_max_list=[];
        con_nomlz_list=[];
    end
    if ~isempty(coneq_list)
        coneq_max_list=max(abs(coneq_list),[],1);
        coneq_nomlz_list=coneq_list./coneq_max_list*10;
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
        % optimal feasiblilty if do not exist feasible point
        object_nonlcon_function_surrogate=@(x) objectNonlconFunctionSurrogate(x,nonlcon_function_surrogate);
        [x_potential,~,exitflag,~]=findMinMSP...
            (object_nonlcon_function_surrogate,variable_number,low_bou,up_bou,[],...
            cheapcon_function,nonlcon_torlance_surrogate);
    end

    % check x_potential if exist in data library
    % if not, updata data libraray
    [x_potential,fval_potential,con_potential,coneq_potential,NFE_p,repeat_index]=dataLibraryUpdataProtect...
        (data_library_name,model_function,x_potential,x_list,low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
    x_list=[x_list;x_potential];
    fval_list=[fval_list;fval_potential];
    if ~isempty(con_list)
        con_list=[con_list;con_potential];
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_potential];
    end

    % normalization fval updata
    if ~isempty(fval_potential)
        fval_potential_nomlz=fval_potential/fval_max*10;
        fval_nomlz_list=[fval_nomlz_list;fval_potential_nomlz];
    end
    if ~isempty(con_potential)
        con_potential_nomlz=(con_potential./con_max_list)*10;
        con_nomlz_list=[con_nomlz_list;con_potential_nomlz];
    end
    if ~isempty(coneq_potential)
        coneq_potential_nomlz=(coneq_potential./coneq_max_list)*10;
        coneq_nomlz_list=[coneq_nomlz_list;coneq_potential_nomlz];
    end

    % when x_potential is exist in data library, x_potential_add will be
    % empty, this times we will use origin point data
    if isempty(x_potential)
        x_potential=x_list(repeat_index,:);
        fval_potential=fval_list(repeat_index,:);
        if ~isempty(con_list)
            con_potential=con_list(repeat_index,:);
        else
            con_potential=[];
        end
        if ~isempty(coneq_list)
            coneq_potential=coneq_list(repeat_index,:);
        else
            coneq_potential=[];
        end
    end

    if DRAW_FIGURE_FLAG && variable_number < 3
        interpVisualize(radialbasis_model_fval,low_bou,up_bou);
        line(x_potential(1),x_potential(2),fval_potential/fval_max*10,'Marker','o','color','r','LineStyle','none')
    end

    % step 6
    % find best result to record
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);

    if INFORMATION_FLAG
        fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
        fprintf('current x:          %s\n',num2str(x_potential));
        fprintf('current value:      %f\n',fval_potential);
        fprintf('current violation:  %s  %s\n',num2str(con_potential),num2str(coneq_potential));
        fprintf('\n');
    end

    result_x_best(iteration,:)=x_best;
    result_fval_best(iteration,:)=fval_best;
    iteration=iteration+1;

    % forced interrupt
    if iteration > iteration_max || NFE >= NFE_max
        done=1;
    end

    % convergence judgment
    if CONVERGENCE_JUDGMENT_FLAG
        if (iteration > 2 && ...
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

    % Interest space sampling
    if ~done

        % step 7
        % using SVM to identify area which is interesting
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
            [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
                (x_list,fval_label,penalty_SVM,low_bou,up_bou);
            if DRAW_FIGURE_FLAG && variable_number < 3
                classifySupportVectorMachineVisualization...
                    (SVM_model,low_bou,up_bou);
            end

            % get data to obtain clustering center
            x_sup_list=[];
            for x_index=1:sample_number_data
                if  SVM_predict_function(x_data_list(x_index,:))==1
                    x_sup_list=[x_sup_list;x_data_list(x_index,:)];
                end
            end

            if isempty(x_sup_list)
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
                    x_center=x_list(filter_index_list(filter_min_index),:);
                else
                    [~,min_fval_index]=min(fval_list(feasible_index_list));
                    x_center=x_list(feasible_index_list(min_fval_index),:);
                end
            end
        else
            % interset sampling
            delta_add=0.1;
            %         fval_thresh=min(fval_list)+(eta-delta_add)*(max(fval_list)-min(fval_list));
            fval_sort=sort(fval_list);
            index=floor(size(fval_list,1)/2);
            fval_thresh=fval_sort(index);
            x_sup_list=[];

            while size(x_sup_list,1) < sample_number_initial/2
                % step 7-1
                % classify exist data
                %             fval_thresh=fval_thresh+delta_add*(max(fval_list)-min(fval_list));
                fval_label=zeros(size(x_list,1),1);
                for x_index=1:size(x_list,1)
                    if fval_list(x_index) <= fval_thresh
                        fval_label(x_index)=1;
                    else
                        fval_label(x_index)=0;
                    end
                end

                % step 7-2
                % get a large number of x point, use SVM to predict x point
                [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
                    (x_list,fval_label,penalty_SVM,low_bou,up_bou);
                if DRAW_FIGURE_FLAG && variable_number < 3
                    classifyVisualization...
                        (SVM_model,low_bou,up_bou);
                end
                for x_index=1:sample_number_data
                    if  SVM_predict_function(x_data_list(x_index,:)')==1
                        x_sup_list=[x_sup_list;x_data_list(x_index,:)];
                    end
                end
            end
        end

        % step 7-3
        % calculate clustering center
        if ~isempty(x_sup_list)
            FC_model=classifyFuzzyClustering...
                (x_sup_list,1,low_bou,up_bou,m);
            x_center=FC_model.center_list;
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
            bou_line=[low_bou_ISR;[low_bou_ISR(1),up_bou_ISR(2)];up_bou_ISR;[up_bou_ISR(1),low_bou_ISR(2)];low_bou_ISR];
            line(bou_line(:,1),bou_line(:,2));
            line(x_potential(1),x_potential(2),'Marker','x')
        end

        % sampling in ISR
        [x_list_exist,~,~,~]=dataLibraryLoad...
            (data_library_name,low_bou_ISR,up_bou_ISR);
        %     [~,x_updata_list,~]=getLatinHypercube...
        %         (sample_number_iteration+size(x_list_exist,1),variable_number,x_list_exist,...
        %         low_bou_ISR,up_bou_ISR,cheapcon_function);
        x_updata_list=(up_bou_ISR-low_bou_ISR).*lhsdesign(sample_number_iteration,variable_number)+low_bou_ISR;

    end

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

    function [x_updata_list,fval_updata_list,con_updata_list,coneq_updata_list,NFE_updata,repeat_index]=dataLibraryUpdataProtect...
            (data_library_name,model_function,x_add_list,...
            x_list,low_bou,up_bou,protect_range)
        % function updata data with same_point_avoid protect
        % return fval
        % all list is x_number x variable_number matrix
        % notice if x_add is exist in library, point will be delete
        %
        variable_number__=size(x_list,2);
        NFE_updata=0;
        x_updata_list=[];fval_updata_list=[];con_updata_list=[];coneq_updata_list=[];repeat_index=[];
        for x_index__=1:size(x_add_list,1)
            x_updata__=x_add_list(x_index__,:);

            % check x_potential if exist in data library
            % if not, updata data libraray
            distance__=sum((abs(x_updata__-x_list)./(up_bou-low_bou)),2);
            [~,min_index__]=min(distance__);
            if distance__(min_index__) < variable_number__*protect_range
                % distance to exist point of point to add is small than protect_range
                repeat_index=[repeat_index;min_index__];
            else
                [fval_updata__,con_updata__,coneq_updata__]=dataLibraryUpdata...
                    (data_library_name,model_function,x_updata__);NFE_updata=NFE_updata+1;
                x_updata_list=[x_updata_list;x_updata__];
                fval_updata_list=[fval_updata_list;fval_updata__];
                con_updata_list=[con_updata_list;con_updata__];
                coneq_updata_list=[coneq_updata_list;coneq_updata__];
            end
        end
    end
end

%% auxiliary function
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
    x=rand(1,variable_number).*(up_bou-low_bou)+low_bou;
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x);
        while sum([~(con < 0);abs(coneq) < 0])
            x=rand(1,variable_number).*(up_bou-low_bou)+low_bou;
            [con,coneq]=cheapcon_function(x);
        end
    end
    population_matrix(population_index,:)=x;
end

% optiaml
ga_option=optimoptions('ga','FunctionTolerance',1e-2,'ConstraintTolerance',1e-2,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',100,'InitialPopulationMatrix',population_matrix,...
    'display','none');
[x_best,fval_best,exitflag,output]=ga...
    (object_function_surrogate,variable_number,[],[],[],[],low_bou',up_bou',constraint_function,ga_option);
fmincon_option=optimoptions('fmincon','FunctionTolerance',1e-6,'ConstraintTolerance',1e-6,...
    'algorithm','sqp',....
    'display','none');
[x_best,fval_best,exitflag,output]=fmincon...
    (object_function_surrogate,x_best,[],[],[],[],low_bou,up_bou,constraint_function,fmincon_option);

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
function FC_model=classifyFuzzyClustering...
    (X,classify_number,low_bou,up_bou,m)
% get fuzzy cluster model
% X is x_number x variable_number matrix
% center_list is classify_number x variable_number matrix
%
if nargin < 4
    up_bou=[];
    if nargin < 3
        low_bou=[];
    end
end
iteration_max=100;
torlance=1e-3;

[x_number,variable_number]=size(X);

% nomalization data
if isempty(low_bou)
    low_bou=min(X,[],1);
else
    low_bou=low_bou(:)';
end
if isempty(up_bou)
    up_bou=max(X,[],1);
else
    up_bou=up_bou(:)';
end
X_nomlz=(X-low_bou)./(up_bou-low_bou);

% if x_number equal 1, clustering cannot done
if x_number==1
    FC_model.X=X;
    FC_model.X_normalize=X_nomlz;
    FC_model.center_list=X;
    FC_model.fval_loss_list=[];
    FC_model.x_class_list=ones(x_number,1);
    return;
end

U=zeros(x_number,classify_number);
center_list=rand(classify_number,variable_number)*0.5;
iteration=0;
done=0;
fval_loss_list=zeros(iteration_max,1);

% get X_center_dis_sq and classify x to each x center
% classify_number x x_number matrix
X_center_dis_sq=zeros(x_number,classify_number);
x_class_list=zeros(x_number,1);
for x_index=1:x_number
    for classify_index=1:classify_number
        temp=(X_nomlz(x_index,:)-center_list(classify_index,:));
        X_center_dis_sq(x_index,classify_index)=temp*temp';
    end
    [~,x_class_list(x_index)]=min(X_center_dis_sq(x_index,:));
end

while ~done
    % updata classify matrix U
    % classify matrix U is x weigth of each center
    for classify_index=1:classify_number
        for x_index=1:x_number
            U(x_index,classify_index)=...
                1/sum((X_center_dis_sq(x_index,classify_index)./X_center_dis_sq(x_index,:)).^(1/(m-1)));
        end
    end

    % updata center_list
    center_list_old=center_list;
    for classify_index=1:classify_number
        center_list(classify_index,:)=...
            sum((U(:,classify_index)).^m.*X_nomlz,1)./...
            sum((U(:,classify_index)).^m,1);
    end

    % updata X_center_dis_sq
    for x_index=1:x_number
        for classify_index=1:classify_number
            temp=(X_nomlz(x_index,:)-center_list(classify_index,:));
            X_center_dis_sq(x_index,classify_index)=temp*temp';
        end
        [~,x_class_list(x_index)]=min(X_center_dis_sq(x_index,:));
    end

    %     plot(center_list(:,1),center_list(:,2));

    % forced interrupt
    if iteration > iteration_max
        done=1;
    end

    % convergence judgment
    if sum(sum(center_list_old-center_list).^2) < torlance
        done=1;
    end

    iteration=iteration+1;
    fval_loss_list(iteration)=sum(sum(U.^m.*X_center_dis_sq));
end
fval_loss_list(iteration+1:end)=[];
center_list=center_list.*(up_bou-low_bou)+low_bou;

FC_model.X=X;
FC_model.X_normalize=X_nomlz;
FC_model.center_list=center_list;
FC_model.fval_loss_list=fval_loss_list;
FC_model.x_class_list=x_class_list;

end

%% SVM
function [predict_function,SVM_model]=classifySupportVectorMachine...
    (X,Class,C,low_bou,up_bou,kernel_function)
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
Y=Class;
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
SVM_model.Class=Class;
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
function [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output]=getRadialBasisModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create radialbasis model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
[predict_function_fval,radialbasis_model_fval]=interpRadialBasisPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    predict_function_con=cell(size(con_list,2),1);
    radialbasis_model_con=struct('X',[],'Y',[],...
        'radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_con=repmat(radialbasis_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        [predict_function_con{con_index},radialbasis_model_con(con_index)]=...
            interpRadialBasisPreModel(x_list,con_list(:,con_index));
    end
else
    predict_function_con=[];
    radialbasis_model_con=[];
end

if ~isempty(coneq_list)
    predict_function_coneq=cell(size(con_list,2),1);
    radialbasis_model_coneq=struct('X',[],'Y',[],...
        'radialbasis_matrix',[],'beta',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],'basis_function',[],...
        'predict_function',[]);
    radialbasis_model_coneq=repmat(radialbasis_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        [predict_function_coneq{coneq_index},radialbasis_model_coneq(coneq_index)]=...
            interpRadialBasisPreModel(x_list,coneq_list(:,coneq_index));
    end
else
    predict_function_coneq=[];
    radialbasis_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,predict_function_fval);
if isempty(predict_function_con) && isempty(predict_function_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,predict_function_con,predict_function_coneq);
end

output.object_function_surrogate=object_function_surrogate;
output.nonlcon_function_surrogate=nonlcon_function_surrogate;

    function fval=objectFunctionSurrogate...
            (X_predict,predict_function_fval)
        % connect all predict favl
        %
        fval=predict_function_fval(X_predict);
    end

    function [con,coneq]=nonlconFunctionSurrogate...
            (X_predict,predict_function_con,predict_function_coneq)
        % connect all predict con and coneq
        %
        if isempty(predict_function_con)
            con=[];
        else
            con=zeros(size(X_predict,1),length(predict_function_con));
            for con_index__=1:length(predict_function_con)
                con(:,con_index__)=....
                    predict_function_con{con_index__}(X_predict);
            end
        end
        if isempty(predict_function_coneq)
            coneq=[];
        else
            coneq=zeros(size(X_predict,1),length(predict_function_coneq));
            for coneq_index__=1:length(predict_function_coneq)
                coneq(:,coneq_index__)=...
                    predict_function_coneq{coneq_index__}(X_predict);
            end
        end
    end
end

function [predict_function,radialbasis_model]=interpRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interp pre model function
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
%
% Copyright 2023 Adel
%
if nargin < 3
    basis_function=[];
end

[x_number,variable_number]=size(X);

% normalize data
aver_X=mean(X);
stdD_X=std(X);
aver_Y=mean(Y);
stdD_Y=std(Y);
index__=find(stdD_X == 0);
if ~isempty(index__),stdD_X(index__)=1;end
index__=find(stdD_Y == 0);
if ~isempty(index__),stdD_Y(index__)=1;end
X_nomlz=(X-aver_X)./stdD_X;
Y_nomlz=(Y-aver_Y)./stdD_Y;

if isempty(basis_function)
    c=(prod(max(X_nomlz)-min(Y_nomlz))/x_number)^(1/variable_number);
    basis_function=@(r) exp(-(r.^2)/c);
end

% initialization distance of all X
X_dis=zeros(x_number,x_number);
for variable_index=1:variable_number
    X_dis=X_dis+(X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end
X_dis=sqrt(X_dis);

[beta,rdibas_matrix]=interpRadialBasis...
    (X_dis,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function=@(X_predict) interpRadialBasisPredictor...
    (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,beta,basis_function);

radialbasis_model.X=X;
radialbasis_model.Y=Y;
radialbasis_model.radialbasis_matrix=rdibas_matrix;
radialbasis_model.beta=beta;

radialbasis_model.aver_X=aver_X;
radialbasis_model.stdD_X=stdD_X;
radialbasis_model.aver_Y=aver_Y;
radialbasis_model.stdD_Y=stdD_Y;
radialbasis_model.basis_function=basis_function;

radialbasis_model.predict_function=predict_function;

% abbreviation:
% num: number, pred: predict, vari: variable
    function [beta,rdibas_matrix]=interpRadialBasis...
            (X_dis,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix=basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix=rdibas_matrix+eye(x_number)*1e-6;

        % solve beta
        beta=rdibas_matrix\Y;
    end

    function [Y_pred]=interpRadialBasisPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,beta,basis_function)
        % radial basis function interpolation predict function
        %
        [x_pred_num,~]=size(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;

        % calculate distance
        X_dis_pred=zeros(x_pred_num,x_num);
        for vari_index=1:vari_num
            X_dis_pred=X_dis_pred+...
                (X_pred_nomlz(:,vari_index)-X_nomlz(:,vari_index)').^2;
        end
        X_dis_pred=sqrt(X_dis_pred);

        % predict variance
        Y_pred=basis_function(X_dis_pred)*beta;

        % normalize data
        Y_pred=Y_pred*stdD_Y+aver_Y;
    end

end

function [ERBF_model_fval,ERBF_model_con,ERBF_model_coneq,output]=getEnsemleRadialBasisModel...
    (x_list,fval_list,con_list,coneq_list)
% base on library_data to create kriging model and function
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
%
[predict_function_fval,ERBF_model_fval]=interpEnsemleRadialBasisPreModel...
    (x_list,fval_list);

if ~isempty(con_list)
    predict_function_con=cell(size(con_list,2),1);
    ERBF_model_con=struct('X',[],'Y',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'basis_function_list',[],'c_list',[],'beta_list',[],...
        'rdibas_matrix_list',[],'inv_rdibas_matrix_list',[],'model_error_list',[],'w',[],...
        'predict_function',[]);
    ERBF_model_con=repmat(ERBF_model_con,[size(con_list,2),1]);
    for con_index=1:size(con_list,2)
        [predict_function_con{con_index},ERBF_model_con(con_index)]=...
            interpolationEnsemleRadialBasisPreModel(x_list,con_list(:,con_index));
    end
else
    predict_function_con=[];
    ERBF_model_con=[];
end

if ~isempty(coneq_list)
    predict_function_coneq=cell(size(con_list,2),1);
    ERBF_model_coneq=struct('X',[],'Y',[],...
        'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
        'basis_function_list',[],'c_list',[],'beta_list',[],...
        'rdibas_matrix_list',[],'inv_rdibas_matrix_list',[],'model_error_list',[],'w',[],...
        'predict_function',[]);
    ERBF_model_coneq=repmat(ERBF_model_coneq,[size(coneq_list,2),1]);
    for coneq_index=1:size(coneq_list,2)
        [predict_function_coneq{con_index},ERBF_model_coneq(con_index)]=...
            interpolationEnsemleRadialBasisPreModel(x_list,coneq_list(:,coneq_index));
    end
else
    predict_function_coneq=[];
    ERBF_model_coneq=[];
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,predict_function_fval);
if isempty(predict_function_con) && isempty(predict_function_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(predict_x) nonlconFunctionSurrogate(predict_x,predict_function_con,predict_function_coneq);
end

output.object_function_surrogate=object_function_surrogate;
output.nonlcon_function_surrogate=nonlcon_function_surrogate;

    function fval=objectFunctionSurrogate...
            (X_predict,predict_function_fval)
        % connect all predict favl
        %
        fval=predict_function_fval(X_predict);
    end

    function [con,coneq]=nonlconFunctionSurrogate...
            (X_predict,predict_function_con,predict_function_coneq)
        % connect all predict con and coneq
        %
        if isempty(predict_function_con)
            con=[];
        else
            con=zeros(size(X_predict,1),length(predict_function_con));
            for con_index__=1:length(predict_function_con)
                con(:,con_index__)=....
                    predict_function_con{con_index__}(X_predict);
            end
        end
        if isempty(predict_function_coneq)
            coneq=[];
        else
            coneq=zeros(size(X_predict,1),length(predict_function_coneq));
            for coneq_index__=1:length(predict_function_coneq)
                coneq(:,coneq_index__)=...
                    predict_function_coneq{coneq_index__}(X_predict);
            end
        end
    end
end

function [predict_function,ensemleradialbasis_model]=interpEnsemleRadialBasisPreModel...
    (X,Y)
% get ensemle radial basis function interpolation model function
% using quadratic programming to calculate weigth of each sub model
% using cubic interpolation optimal to decrese time use
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
%
% reference: [1] SHI R, LIU L, LONG T, et al. An efficient ensemble of
% radial basis functions method based on quadratic programming [J].
% Engineering Optimization, 2016, 48(1202 - 25.
%
% Copyright 2023 Adel
%
[x_number,variable_number]=size(X);

% normalize data
aver_X=mean(X);
stdD_X=std(X);
aver_Y=mean(Y);
stdD_Y=std(Y);
index__=find(stdD_X == 0);
if ~isempty(index__),stdD_X(index__)=1;end
index__=find(stdD_Y == 0);
if ~isempty(index__),stdD_Y(index__)=1;end
X_nomlz=(X-aver_X)./stdD_X;
Y_nomlz=(Y-aver_Y)./stdD_Y;

c_initial=(prod(max(X_nomlz)-min(Y_nomlz))/x_number)^(1/variable_number);

% initialization distance of all X
X_dis=zeros(x_number,x_number);
for variable_index=1:variable_number
    X_dis=X_dis+(X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end
X_dis=sqrt(X_dis);

% linear kernal function
basis_function_linear=@(r,c) r+c;
dRM_dc_function=@(x_number,X_dis,rdibas_matrix,c) ones(x_number,x_number);
object_function=@(c) objectFunctionRadiabasis....
    (X_dis,Y_nomlz,x_number,basis_function_linear,c,dRM_dc_function);

[linear_c_backward,fval_backward,~]=optimalCubicInterp...
    (object_function,-1e2,-1e2,1e2,1e-3);
[linear_c_forward,fval_forward,~]=optimalCubicInterp...
    (object_function,1e2,-1e2,1e2,1e-3);
if fval_forward < fval_backward
    c_linear=linear_c_forward;
else
    c_linear=linear_c_backward;
end

% gauss kernal function
basis_function_gauss=@(r,c) exp(-c*r.^2);
dRM_dc_function=@(x_number,X_dis,rdibas_matrix,c) -X_dis.^2.*rdibas_matrix;
object_function=@(c) objectFunctionRadiabasis....
    (X_dis,Y_nomlz,x_number,basis_function_gauss,c,dRM_dc_function);
[c_gauss,~,~,~]=optimalCubicInterp...
    (object_function,c_initial,1e-2,1e2,1e-3);

% spline kernal function
basis_function_spline=@(r,c) r.^2.*log(r.^2*c+1e-3);
dRM_dc_function=@(x_number,X_dis,rdibas_matrix,c) X_dis.^4./(X_dis.^2*c+1e-3);
object_function=@(c) objectFunctionRadiabasis....
    (X_dis,Y_nomlz,x_number,basis_function_spline,c,dRM_dc_function);
[c_spline,~,~,~]=optimalCubicInterp...
    (object_function,c_initial,1e-2,1e2,1e-3);

% triple kernal function
basis_function_triple=@(r,c) (r+c).^3;
dRM_dc_function=@(x_number,X_dis,rdibas_matrix,c) 3*(X_dis+c).^3;
object_function=@(c) objectFunctionRadiabasis....
    (X_dis,Y_nomlz,x_number,basis_function_triple,c,dRM_dc_function);
[c_triple,~,~,~]=optimalCubicInterp...
    (object_function,c_initial,1e-2,1e2,1e-3);

% multiquadric kernal function
basis_function_multiquadric=@(r,c) sqrt(r+c);
dRM_dc_function=@(x_number,X_dis,rdibas_matrix,c) 0.5./rdibas_matrix;
object_function=@(c) objectFunctionRadiabasis....
    (X_dis,Y_nomlz,x_number,basis_function_multiquadric,c,dRM_dc_function);
[c_binomial,~,~,~]=optimalCubicInterp...
    (object_function,c_initial,1e-2,1e2,1e-3);

% inverse multiquadric kernal function
basis_function_inverse_multiquadric=@(r,c) 1./sqrt(r+c);
dRM_dc_function=@(x_number,X_dis,rdibas_matrix,c) -0.5*rdibas_matrix.^3;
object_function=@(c) objectFunctionRadiabasis....
    (X_dis,Y_nomlz,x_number,basis_function_inverse_multiquadric,c,dRM_dc_function);

% c_initial=1;
% drawFunction(object_function,1e-1,10);

[c_inverse_binomial,~,~,~]=optimalCubicInterp...
    (object_function,c_initial,1e-2,1e2,1e-3);

% generate total model
basis_function_list={
    @(r) basis_function_linear(r,c_linear);
    @(r) basis_function_gauss(r,c_gauss);
    @(r) basis_function_spline(r,c_spline);
    @(r) basis_function_triple(r,c_triple);
    @(r) basis_function_multiquadric(r,c_binomial);
    @(r) basis_function_inverse_multiquadric(r,c_inverse_binomial);};
c_list=[c_linear;c_gauss;c_spline;c_triple;c_binomial;c_inverse_binomial];

model_number=size(basis_function_list,1);
beta_list=zeros(x_number,model_number);
rdibas_matrix_list=zeros(x_number,x_number,model_number);
inv_rdibas_matrix_list=zeros(x_number,x_number,model_number);

% calculate model matrix and error
model_error_list=zeros(model_number,x_number);
for model_index=1:model_number
    basis_function=basis_function_list{model_index};
    [beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
        (X_dis,Y_nomlz,basis_function,x_number);
    beta_list(:,model_index)=beta;
    rdibas_matrix_list(:,:,model_index)=rdibas_matrix;
    inv_rdibas_matrix_list(:,:,model_index)=inv_rdibas_matrix;

    model_error_list(model_index,:)=(beta./...
        diag(inv_rdibas_matrix))';
end

% calculate weight of each model
C=model_error_list*model_error_list';
eta=trace(C)/x_number;
I_model=eye(model_number);
one_model=ones(model_number,1);
w=(C+eta*I_model)\one_model/...
    (one_model'*((C+eta*I_model)\one_model));
while min(w) < -0.05
    % minimal weight cannot less than zero too much
    eta=eta*10;
    w=(C+eta*I_model)\one_model/...
        (one_model'*((C+eta*I_model)\one_model));
end

% initialization predict function
predict_function=@(X_predict) interpEnsemleRadialBasisPredictor...
    (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,model_number,beta_list,basis_function_list,w);

ensemleradialbasis_model.X=X;
ensemleradialbasis_model.Y=Y;
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

% abbreviation:
% num: number, pred: predict, vari: variable
    function [fval,gradient]=objectFunctionRadiabasis....
            (X_dis,Y,x_number,basis_function,c,dRM_dc_function)
        % MSE_CV function, simple approximation to RMS
        % basis_function input is c and x_sq
        %
        basis_function=@(r) basis_function(r,c);
        [beta__,rdibas_matrix__,inv_rdibas_matrix__]=interpRadialBasis...
            (X_dis,Y,basis_function,x_number);
        U=beta__./diag(inv_rdibas_matrix__);
        fval=sum(U.^2);

        % calculate gradient
        if nargout > 1
            inv_rdibas_matrix_gradient=-inv_rdibas_matrix__*...
                dRM_dc_function...
                (x_number,X_dis,rdibas_matrix__,c)*inv_rdibas_matrix__;
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

    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
            (X_dis,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix=basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix=rdibas_matrix+eye(x_number)*1e-6;

        % solve beta
        inv_rdibas_matrix=rdibas_matrix\eye(x_number);
        beta=inv_rdibas_matrix*Y;
    end

    function [x_best,favl_best,NFE,output]=optimalCubicInterp...
            (object_function,x_initial,low_bou,up_bou,torlance,iteration_max)
        % cubic interp optimization, should provide fval and gradient
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

        % main loop for cubic interp
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
            if abs(x-x_old) < 1e-5
                quit_flag=1;
            end
            if iteration >= iteration_max
                quit_flag=1;
            end
        end
    end

    function Y_pred=interpEnsemleRadialBasisPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,model_number,beta_list,basis_function_list,w)
        % ensemle radial basis function interpolation predict function
        %
        [x_pred_num,~]=size(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;

        % calculate distance
        X_dis_pred=zeros(x_pred_num,x_num);
        for vari_index=1:vari_num
            X_dis_pred=X_dis_pred+...
                (X_pred_nomlz(:,vari_index)-X_nomlz(:,vari_index)').^2;
        end
        X_dis_pred=sqrt(X_dis_pred);

        % calculate each sub model predict fval and get predict_y
        y_pred_nomlz_list=zeros(x_pred_num,model_number);
        for model_index__=1:model_number
            basis_function__=basis_function_list{model_index__};
            beta__=beta_list(:,model_index__);
            y_pred_nomlz_list(:,model_index__)=basis_function__(X_dis_pred)*beta__;
        end
        Y_pred_nomlz=y_pred_nomlz_list*w;

        % normalize data
        Y_pred=Y_pred_nomlz*stdD_Y+aver_Y;
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
    fval_list=[fval_list;fval(:)'];
    con_list=[con_list;con(:)'];
    coneq_list=[coneq_list;coneq(:)'];

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
                    distance_min__=distance_temp__;
                end
                if distance_temp__ < search_range__
                    search_range__=distance_temp__;
                end
                x_next_index__=x_next_index__+1;
            end
        end
        distance_min__=sqrt(distance_min__);
    end
end
