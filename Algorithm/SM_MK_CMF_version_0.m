clc;
clear;
close all hidden;

benchmark_function=BenchmarkFunction();

% variable_number=2;
% object_function=@(x) functionGPObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-2;-2];
% up_bou=[2;2];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

% variable_number=2;
% object_function=@(x) functionBRObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-5;10];
% up_bou=[0;15];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

% variable_number=2;
% object_function=@(x) functionSCObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-2;-2];
% up_bou=[2;2];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

% variable_number=2;
% object_function=@(x) functionRSObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-1;-1];
% up_bou=[1;1];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

% variable_number=2;
% object_function=@(x) functionPKObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-3;-3];
% up_bou=[3;3];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

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

variable_number=2;
object_function=@functionG06Object;
A=[];
B=[];
Aeq=[];
Beq=[];
low_bou=[13;0];
up_bou=[100;100];
cheapcon_function=[];
nonlcon_function=@functionG06Nonlcon;
model_function=[];

% variable_number=4;
% object_function=@(x) functionPVD4Object(x);
% object_function_low=@(x) functionPVD4ObjectLow(x);
% A=[-1,0,0.0193,0;
%     0,-1,0.00954,0;];
% B=[0;0];
% Aeq=[];
% Beq=[];
% low_bou=[0;0;0;0];
% up_bou=[1;1;50;240];
% nonlcon_function=@(x) functionPVD4Nonlcon(x);
% cheapcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq);
% model_function=[];

data_library_name='optimal_SM_MK_CMF_result.txt';
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

function [con,coneq]=cheapconFunction(x,A,B,Aeq,Beq,cheapcon_function)
% convert A, B, Aeq, Beq to total cheapcon function
%
if nargin < 6
    cheapcon_function=[];
    if nargin < 5
        Beq=[];
        if nargin < 4
            Aeq=[];
            if nargin < 3
                B=[];
                if nargin < 2
                    A=[];
                end
            end
        end
    end
end
x=x(:);
con=[];
coneq=[];
if ~isempty(A)
    if isempty(B)
        con=[con;A*x];
    else
        con=[con;A*x-B];
    end
end
if ~isempty(Aeq)
    if isempty(Beq)
        coneq=[coneq;Aeq*x];
    else
        coneq=[coneq;Aeq*x-Beq];
    end
end
if ~isempty(cheapcon_function)
    [lincon,linconeq]=cheapcon_function(x);
    con=[con;lincon];
    coneq=[coneq;linconeq];
end
end

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

data_library_name='optimal_SM_MK_CMF_result';
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

function CGPMF_model=classifyGaussProcessMultiFidelity...
    (X_L,Y_L,X_H,Y_H,low_bou,up_bou)
% generate gauss classifier model
% X is x_number x variable_number matirx,Y is x_number x 1 matrix
% low_bou, up_bou is variable_number x 1 matrix
% support multi fidelity
%
if nargin < 4
    up_bou=[];
    if nargin < 3
        low_bou=[];
    end
end

[x_L_number,variable_number]=size(X_L);
[x_H_number,~]=size(X_H);
x_number=x_L_number+x_H_number;

% poster function list, x_number f, variable_number len, eta
poster_function_list=cell(x_number+2*(variable_number+1)+1,1);
for x_index=1:x_number
    poster_function_list{x_index}=@(x) logNormalFunction(x,0,1);
end
base=x_number;
for variable_index=1:variable_number
    poster_function_list{base+variable_index}=@(x) logGammaFunction(x,2,2);
end
poster_function_list{base+variable_number+1}=@(x) logHalfNormalFunction(x,5);
base=x_number+variable_number+1;
for variable_index=1:variable_number
    poster_function_list{base+variable_index}=@(x) logGammaFunction(x,2,2);
end
poster_function_list{base+variable_number+1}=@(x) logHalfNormalFunction(x,5);
poster_function_list{end}=@(x) logNormalFunction(x,0,10);

% normalization data
if isempty(up_bou)
    up_bou=max([X_L;X_H])';
end
if isempty(low_bou)
    low_bou=min([X_L;X_H])';
end
X_nomlz=(X_L-low_bou')./(up_bou'-low_bou');
X_H_nomlz=(X_H-low_bou')./(up_bou'-low_bou');

% total list
X=[X_L;X_H];
X_nomlz=[X_nomlz;X_H_nomlz];
Y=[Y_L;Y_H];

% initializate square of X inner distance x__x, x__x_h, x_h__x_h sq total
X_dis_sq=zeros(x_number,x_number,variable_number);
for index_i=1:x_L_number
    % x__x
    for index_j=1:index_i-1
        X_dis_sq(index_i,index_j,:)=X_dis_sq(index_j,index_i,:);
    end
    X_dis_sq(index_i,index_i,:)=0;
    for index_j=index_i+1:x_L_number
        X_dis_sq(index_i,index_j,:)=(X_nomlz(index_i,:)-X_nomlz(index_j,:)).^2;
    end
    % x__x_h
    for index_j=x_L_number+1:x_number
        X_dis_sq(index_i,index_j,:)=(X_nomlz(index_i,:)-X_nomlz(index_j,:)).^2;
    end
end
for index_i=x_L_number+1:x_number
    % x_h__x and x_h__x_h
    for index_j=1:index_i-1
        X_dis_sq(index_i,index_j,:)=X_dis_sq(index_j,index_i,:);
    end
    X_dis_sq(index_i,index_i,:)=0;
    % x_h__x_h
    for index_j=x_L_number+1:x_number
        X_dis_sq(index_i,index_j,:)=(X_nomlz(index_i,:)-X_nomlz(index_j,:)).^2;
    end
end

% likelihood function
log_likelihood_function=@(hyperparameter) logLikelihoodFunction...
    (hyperparameter,X_dis_sq,Y,...
            x_L_number,x_H_number,x_number,variable_number);
% posterior function
log_posterior_function=@(hyperparameter) logPosteriorFunction...
    (log_likelihood_function,hyperparameter,poster_function_list);

% get hyperparameter
hyperparameter_initial=ones(x_number+2*(variable_number+1)+1,1);

% [fval,gradient]=log_likelihood_function(hyperparameter_initial)
% [gradient_differ,~]=differ...
%     (log_likelihood_function,hyperparameter_initial,fval)

low_bou_hyper=-realmax*ones(x_number+2*(variable_number+1)+1,1);
low_bou_hyper(x_number+(variable_number+1))=1e-6;
low_bou_hyper(x_number+2*(variable_number+1))=1e-6;
up_bou_hyper=realmax*ones(x_number+2*(variable_number+1)+1,1);

% optimal method
object_function=@(x) objectFunction(log_posterior_function,x);
fmincon_option=optimoptions('fmincon','Display','none','Algorithm','sqp','SpecifyObjectiveGradient',true);
[hyperparameter,fval,~,~]=fmincon...
    (object_function,hyperparameter_initial,[],[],[],[],low_bou_hyper,up_bou_hyper,[],fmincon_option);

% % sampling method to get hyperparamter
% sample_number=400;
% iteration_adapt=400;
% accept_probability=0.65;
% [sample_list,output]=samplerNUTS...
%     (log_posterior_function,hyperparameter_initial,sample_number,iteration_adapt,accept_probability,...
%     low_bou_hyper,up_bou_hyper);
% hyperparameter=(sum(sample_list(:,:),1)/(sample_number))';
% count_number=100;
% [~,~,~,~]=drawSample...
%     (sample_list,count_number);
% CGPMF_model.sample_list=sample_list;

disp(['hyperparameter: ']);disp(hyperparameter(x_number+1:end));

v=hyperparameter(1:x_number);
base=x_number;
len_L=hyperparameter((base+1):(base+variable_number));
eta_L=hyperparameter(base+variable_number+1);
base=x_number+variable_number+1;
len_H=hyperparameter((base+1):(base+variable_number));
eta_H=hyperparameter(base+variable_number+1);
rou=hyperparameter(end);
% obtain convariance and inv_covariance_X
[covariance_X,inv_covariance_X,...
    ~,~,~,...
    ~,~]=getCovariance...
    (len_L,eta_L,len_H,eta_H,rou,...
    X_dis_sq,x_L_number,x_H_number,x_number,variable_number);
fvpre=chol(covariance_X,'lower')*v;

% generate predict function
CGMFM_predict_function=@(x) classifyGaussPredictor...
    (x,X_nomlz,fvpre,len_L,eta_L,len_H,eta_H,rou,inv_covariance_X,...
    low_bou,up_bou,x_L_number,x_H_number,x_number,variable_number);

% output model
CGPMF_model.X=X_L;
CGPMF_model.Y=Y_L;
CGPMF_model.X_H=X_H;
CGPMF_model.Y_H=Y_H;
CGPMF_model.X_total=X;
CGPMF_model.Y_total=Y;
CGPMF_model.X_nomlz=X_nomlz;
CGPMF_model.X_H_nomlz=X_H_nomlz;
CGPMF_model.low_bou=low_bou;
CGPMF_model.up_bou=up_bou;
CGPMF_model.inv_covariance_X=inv_covariance_X;
CGPMF_model.CGMFM_predict_function=CGMFM_predict_function;

    function [fval,gradient]=logLikelihoodFunction...
            (hyperparameter,X_dis_sq,Y,...
            x_L_number,x_H_number,x_number,variable_number)
        % support vector machine maximum object function
        % hyperparameter is x_number v,
        % variable_number len_L, eta_L, variable_number len_H, eta_H rou
        % only calculate half of matrix, due to symmetry of matrix
        %
        v__=hyperparameter(1:x_number);
        base__=x_number;
        len_L__=hyperparameter((base__+1):(base__+variable_number));
        eta_L__=hyperparameter(base__+variable_number+1);
        base__=x_number+variable_number+1;
        len_H__=hyperparameter((base__+1):(base__+variable_number));
        eta_H__=hyperparameter(base__+variable_number+1);
        rou__=hyperparameter(end);
        
        [covariance_X__,inv_covariance_X__,...
            exp_dis__,rou_exp_dis__,eta_rou_exp_dis__,...
            exp_H_dis__,eta_H_exp_H_dis__]=getCovariance...
            (len_L__,eta_L__,len_H__,eta_H__,rou__,...
            X_dis_sq,x_L_number,x_H_number,x_number,variable_number);
        
        L=chol(covariance_X__,'lower');
        inv_L=inv(L);
        fvpre__=L*v__;
        exp__fvpre=exp(-fvpre__);
        p__=(1-2e-6)./(1+exp__fvpre)+1e-6;
        log_p__=log(p__);
        log_1_p__=log(1-p__);
        
        fval=sum(Y.*log_p__+(1-Y).*log_1_p__);
        
        % get gradient
        dfval_dfvpre=(Y./p__-(1-Y)./(1-p__)).*exp__fvpre./((1+exp__fvpre).^2);
        gradient=zeros(x_number+2*(variable_number+1)+1,1);
        % fvpre
        for x_index__=1:x_number
            gradient(x_index__)=sum(dfval_dfvpre.*L(:,x_index__));
        end
        base__=x_number;
        % len_L
        for variable_index__=1:variable_number
            temp=len_L__(variable_index__)^3;
            % cov_gradient is dcov_dlen
            cov_gradient=X_dis_sq(:,:,variable_index__).*eta_rou_exp_dis__/temp;
            dcfp_dlen_L=(L*(semidiagonal(x_number).*(inv_L*cov_gradient*inv_L')))*v__;
            gradient(base__+variable_index__)=sum(dfval_dfvpre.*dcfp_dlen_L);
        end
        % eta_L
        % cov_gradient=rou_exp_dis__ is dcov_deta
        dcfp_deta_L=(L*(semidiagonal(x_number).*(inv_L*rou_exp_dis__*inv_L')))*v__;
        gradient(base__+variable_number+1)=sum(dfval_dfvpre.*dcfp_deta_L);
        
        % rou
        exp_dis__(1:x_L_number,1:x_L_number)=0;
        cov_gradient=exp_dis__;
        cov_gradient(1:x_L_number,x_L_number+1:x_number)=eta_L__*...
            cov_gradient(1:x_L_number,x_L_number+1:x_number);
        cov_gradient(x_L_number+1:x_number,1:x_L_number)=...
            cov_gradient(1:x_L_number,x_L_number+1:x_number)';
        cov_gradient(x_L_number+1:x_number,x_L_number+1:x_number)=2*rou__*eta_L__*...
            cov_gradient(x_L_number+1:x_number,x_L_number+1:x_number);
        dcfp_drou=(L*(semidiagonal(x_number).*(inv_L*cov_gradient*inv_L')))*v__;
        gradient(end)=sum(dfval_dfvpre.*dcfp_drou);
        
        % len_H
        base__=x_number+variable_number+1;
        X_dis_sq(1:x_L_number,1:x_L_number)=0;
        X_dis_sq(1:x_L_number,x_L_number+1:x_number)=0;
        X_dis_sq(x_L_number+1:x_number,1:x_L_number)=0;
        
        for variable_index__=1:variable_number
            temp=len_H__(variable_index__)^3;
            % cov_gradient is dcov_dlen
            cov_gradient=X_dis_sq(:,:,variable_index__).*eta_H_exp_H_dis__/temp;
            dcfp_dlen_H=(L*(semidiagonal(x_number).*(inv_L*cov_gradient*inv_L')))*v__;
            gradient(base__+variable_index__)=sum(dfval_dfvpre.*dcfp_dlen_H);
        end
        % eta_H
        % cov_gradient=exp_H_dis__ is dcov_deta
        dcfp_deta_H=(L*(semidiagonal(x_number).*(inv_L*exp_H_dis__*inv_L')))*v__;
        gradient(base__+variable_number+1)=sum(dfval_dfvpre.*dcfp_deta_H);
        
        function matrix=semidiagonal(n)
            matrix=ones(n);
            for rank_index=1:n
                matrix(rank_index,rank_index)=0.5;
                for colume_index=rank_index+1:n
                    matrix(rank_index,colume_index)=0;
                end
            end
        end
    end
    function [fval,gradient]=logPosteriorFunction...
            (log_likelihood_function,hyperparameter,poster_function_list)
        % add poster disturbution
        %
        [fval,gradient]=log_likelihood_function(hyperparameter);
        hyper_number=size(hyperparameter,1);
        for hyper_index=1:hyper_number
            poster_function=poster_function_list{hyper_index};
            [fval_hyper,gradient_hyper]=poster_function(hyperparameter(hyper_index));
            fval=fval+fval_hyper;
            gradient(hyper_index)=gradient(hyper_index)+gradient_hyper;
        end
    end
    function [covariance_X__,inv_covariance_X__,...
            exp_dis__,rou_exp_dis__,eta_rou_exp_dis__,...
            exp_H_dis__,eta_H_exp_H_dis__]=getCovariance...
            (len_L__,eta_L__,len_H__,eta_H__,rou__,...
            X_dis_sq,x_L_number,x_H_number,x_number,variable_number)
        % obtain covariance of x
        rou_sq__=rou__*rou__;
        
        % exp of x__x, x__x_h, x_h__x_h with theta
        exp_dis__=zeros(x_number);
        for rank_index__=1:x_number
            % symmetry
            for colume_index__=1:rank_index__-1
                exp_dis__(rank_index__,colume_index__)=exp_dis__(colume_index__,rank_index__);
            end
            
            % diagonal
            exp_dis__(rank_index__,rank_index__)=1+1e-3;
            
            % initial
            for colume_index__=rank_index__+1:x_number
                temp=0;
                for variable_index__=1:variable_number
                    temp=temp+X_dis_sq(rank_index__,colume_index__,variable_index__)/...
                        (2*len_L__(variable_index__)^2);
                end
                exp_dis__(rank_index__,colume_index__)=exp(-temp);
            end
        end
        
        % only exp of x_h__x_h with theta_H
        exp_H_dis__=zeros(x_number);
        for rank_index__=x_L_number+1:x_number
            % symmetry
            for colume_index__=1:rank_index__-1
                exp_H_dis__(rank_index__,colume_index__)=exp_H_dis__(colume_index__,rank_index__);
            end
            
            % diagonal
            exp_H_dis__(rank_index__,rank_index__)=1+1e-3;
            
            % initial
            for colume_index__=rank_index__+1:x_number
                temp=0;
                for variable_index__=1:variable_number
                    temp=temp+X_dis_sq(rank_index__,colume_index__,variable_index__)/...
                        (2*len_H__(variable_index__)^2);
                end
                exp_H_dis__(rank_index__,colume_index__)=exp(-temp);
            end
        end
        
        % add rou
        rou_exp_dis__=exp_dis__;
        for rank_index__=x_L_number+1:x_number
            % x__x_h and symmetry
            for colume_index__=1:x_L_number
                rou_exp_dis__(rank_index__,colume_index__)=rou_exp_dis__(rank_index__,colume_index__)*rou__;
                rou_exp_dis__(colume_index__,rank_index__)=rou_exp_dis__(rank_index__,colume_index__);
            end
            
            % symmetry
            for colume_index__=x_L_number+1:rank_index__-1
                rou_exp_dis__(rank_index__,colume_index__)=rou_exp_dis__(colume_index__,rank_index__);
            end
            % x_h__x_h
            for colume_index__=rank_index__:x_number
                rou_exp_dis__(rank_index__,colume_index__)=rou_exp_dis__(rank_index__,colume_index__)*rou_sq__;
            end
        end
        
        eta_rou_exp_dis__=rou_exp_dis__*eta_L__;
        eta_H_exp_H_dis__=exp_H_dis__*eta_H__;
        
        % convariance of total data
        covariance_X__=eta_rou_exp_dis__+eta_H_exp_H_dis__;
        
        inv_covariance_X__=inv(covariance_X__);
    end
    function [fval,gradient]=objectFunction(initial_function,x)
        [fval,gradient]=initial_function(x);
        fval=-fval;
        gradient=-gradient;
    end
    function [class,mu_pre,var_pre]=classifyGaussPredictor...
            (x,X_nomlz,fvpre,len_L,eta_L,len_H,eta_H,rou,inv_covariance_X,...
            low_bou,up_bou,x_L_number,x_H_number,x_number,variable_number)
        % predict value of x
        % x input is colume vector
        % X_nomlz include X_L and X_H
        %
        x_nomlz=(x-low_bou)./(up_bou-low_bou);
        x_cov__=zeros(x_number,1);
        % cov of x and X_L
        for x_index__=1:x_L_number
            temp=sum((x_nomlz'-X_nomlz(x_index__,:)).^2./(len_L').^2/2);
            x_cov__(x_index__)=rou*eta_L*exp(-temp);
        end
        % cov of x and X_H
        for x_index__=x_L_number+1:x_number
            temp=sum((x_nomlz'-X_nomlz(x_index__,:)).^2./(len_L').^2/2);
            temp_H=sum((x_nomlz'-X_nomlz(x_index__,:)).^2./(len_H').^2/2);
            x_cov__(x_index__)=rou*rou*eta_L*exp(-temp)+eta_H*exp(-temp_H);
        end
        
        % get mu_pre
        mu_pre=x_cov__'*inv_covariance_X*fvpre;
        if 1/(1+exp(-mu_pre)) > 0.5
            class=1;
        else
            class=0;
        end
        var_pre=rou*rou*eta_L+eta_H-x_cov__'*inv_covariance_X*x_cov__;
    end

end
function classifyGaussProcessMultiFidelityVisualization...
    (CGPMF_model,low_bou,up_bou,figure_handle)
% Visualization SVM_model
%
if nargin < 1
    error('classifySupportVectorMachineVisualization: not enough input');
end
X=CGPMF_model.X_H;
Y=CGPMF_model.Y_H;

CGMFM_predict_function=CGPMF_model.CGMFM_predict_function;

if nargin < 3
    up_bou=max(X)';
    if nargin < 2
        low_bou=min(X)';
    end
end

if nargin < 4
    figure_handle=figure(200);
end
axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end

% check dimension
[~,variable_number]=size(X);
if variable_number ~= 2
    error('classifySupportVectorMachineVisualization: dimension do not equal 2');
end

index=find(Y>0);
X_positive=X(index,:);
index=find(Y<0);
X_negative=X(index,:);

% draw zero value line
grid_number=100;
d_bou=(up_bou-low_bou)/grid_number;
[X_draw,Y_draw]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
class=zeros(grid_number+1);
fval=zeros(grid_number+1);
for x_index=1:grid_number+1
    for y_index=1:grid_number+1
        predict_x=([x_index;y_index]-1).*d_bou+low_bou;
        [class(y_index,x_index),fval(y_index,x_index)]=...
            CGMFM_predict_function(predict_x);
    end
end
contour(axes_handle,X_draw,Y_draw,class);

% draw point
line(axes_handle,X_positive(:,1),X_positive(:,2),'LineStyle','none','Marker','o','Color','b');
line(axes_handle,X_negative(:,1),X_negative(:,2),'LineStyle','none','Marker','o','Color','r');

% figure(2)
% surf(X,Y,fval_sum);
end

function [fval,gradient]=logHalfNormalFunction(x,sigma)
if x < 0
    fval=-realmax;
    gradient=1e6;
else
    sigma_sq=sigma*sigma;
    fval=-0.5*log(pi/2)-log(sigma)-x*x/2/sigma_sq;
    gradient=-x/sigma_sq;
end
end
function [fval,gradient]=logNormalFunction(x,mu,sigma)
sigma_sq=sigma*sigma;
x_mu=x-mu;
fval=-0.5*log(2*pi)-log(sigma)-x_mu*x_mu/2/sigma_sq;
gradient=-x_mu/sigma_sq;
end
function [fval,gradient]=logGammaFunction(x,alpha,beta)
if x < 0
    fval=-realmax;
    gradient=1e6;
else
    fval=-log(beta^alpha)-log(gamma(alpha))+(alpha-1)*log(x)-x/beta;
    gradient=(alpha-1)/x-1/beta;
end
end

function [MK_model_fval,MK_model_con,MK_model_coneq,output]=getMKModel...
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

function [predict_function,MK_model]=interpMultiKriging...
    (XHF, YHF, XLF, YLF)
% construct Multi-Level Kriging 
% XHF, YHF are x_HF_number x variable_number matrix
% XLF, YLF are x_LF_number x variable_number matrix
% aver_X,stdD_X is 1 x x_HF_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
%
% input:
% XHF, YHF, XLF, YLF
%
% output:
% kriging model (predict_function,
% XHF, YHF, XLF_list, YLF_list, base_function_list)
%
% Copyright 2022.12 Adel
%
[x_HF_number,variable_number]=size(XHF);

% construct low fidelity
[predict_function_LF,kriging_model_LF]=interpKrigingPreModel...
    (XLF,YLF);

% evaluate error in high fidelity point
LF_predict_list=zeros(x_HF_number,1); % 
for x_index=1:x_HF_number
    LF_predict_list(x_index,1)=...
        kriging_model_LF.predict_function(XHF(x_index,:));
end

% construct bias kriging model
Y_bias=YHF-LF_predict_list;
[predict_function_bias,kriging_model_bias]=interpKrigingPreModel...
    (XHF,Y_bias);

% initialization predict function
predict_function=@(predict_x) interpMultiKrigingPredictor...
    (predict_x,predict_function_LF,predict_function_bias);

MK_model.XHF=XHF;
MK_model.X=XHF;
MK_model.YHF=YHF;
MK_model.Y=YHF;
MK_model.kriging_model_LF=kriging_model_LF;
MK_model.kriging_model_bias=kriging_model_bias;
MK_model.predict_function=predict_function;

    function predict_fval=interpMultiKrigingPredictor...
            (predict_x,predict_function_LF,predict_function_bias)
        % kriging interpolation predict function
        %
        predict_fval=predict_function_LF(predict_x)+...
            predict_function_bias(predict_x);
    end
end

function [predict_function,kriging_model]=interpKrigingPreModel...
    (X,Y,theta)
% version 4, nomalization method is grassian
% prepare model, optimal theta and calculation parameter
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
%
% input initial data X, Y, which are real data
%
% output is a kriging model, include predict_function...
% X, Y, base_function_list
%
% Copyright 2022.10 Adel
%
[x_number,variable_number]=size(X);
if nargin < 3
    theta=ones(1,variable_number);
end

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__=find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__)=1; end
index__=find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__)=1; end
X_nomlz=(X-repmat(aver_X,x_number,1))./repmat(stdD_X,x_number,1);
Y_nomlz=(Y-repmat(aver_Y,x_number,1))./repmat(stdD_Y,x_number,1);

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
fmincon_option=optimoptions(@fmincon,'Display','iter-detailed',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10);
low_bou_kriging=1e-1*ones(variable_number,1);
up_bou_kriging=20*ones(variable_number,1);
object_function_kriging=@(theta) objectFunctionKriging...
    (X_dis_sq,X_nomlz,Y_nomlz,x_number,variable_number,theta);

theta=fmincon...
    (object_function_kriging,theta,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% get parameter
[covariance,inv_covariance,fval_reg,beta,sigma_sq]=interpKriging...
    (X_dis_sq,X_nomlz,Y_nomlz,x_number,variable_number,theta);
gama=inv_covariance*(Y_nomlz-fval_reg*beta);
FTRF=fval_reg'*inv_covariance*fval_reg;

% initialization predict function
predict_function=@(predict_x) interpKrigingPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    theta,beta,gama,sigma_sq,...
    inv_covariance,fval_reg,FTRF,predict_x);

kriging_model.X=X;
kriging_model.Y=Y;
kriging_model.X_normalize=X_nomlz;
kriging_model.Y_normalize=Y_nomlz;
kriging_model.fval_regression=fval_reg;
kriging_model.covariance=covariance;
kriging_model.inv_covariance=inv_covariance;

kriging_model.theta=theta;
kriging_model.beta=beta;
kriging_model.gama=gama;
kriging_model.sigma_sq=sigma_sq;
kriging_model.aver_X=aver_X;
kriging_model.stdD_X=stdD_X;
kriging_model.aver_Y=aver_Y;
kriging_model.stdD_Y=stdD_Y;

kriging_model.predict_function=predict_function;

    function fval=objectFunctionKriging...
            (X_dis_sq,X,Y,x_number,variable_number,theta)
        % function to minimize sigma_sq
        %
        cov__=zeros(x_number,x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                cov__(rank_index__,colume_index__)=...
                    cov__(colume_index__,rank_index__);
            end
            cov__(rank_index__,rank_index__)=1+1e-6; % stabilize
            for colume_index__=rank_index__+1:x_number
                temp__=X_dis_sq(rank_index__,colume_index__,:);
                cov__(rank_index__,colume_index__)=...
                    exp(-temp__(:)'*theta(:));
            end
        end
        
        % F is base funcion fval of data_point_x
        fval_reg__=[ones(x_number,1),X];
        
        % coefficient calculation
        inv_cov__=inv(cov__);
        beta__=(fval_reg__'*inv_cov__*fval_reg__)\fval_reg__'*inv_cov__*Y;
        Y_fbeta=Y-fval_reg__*beta__;
        fval=(Y_fbeta'*inv_cov__*Y_fbeta)/x_number;
        
    end
    function [cov,inv_cov,fval_reg,beta,sigma_sq]=interpKriging...
            (X_dis_sq,X,Y,x_number,variable_number,theta)
        % total riging interpolation function
        %
        % input X, Y as initial data, theta and base function
        %
        % output covariance, inv of covariance,...
        %
        % Copyright 2022 Adel
        %
        cov=zeros(x_number,x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                cov(rank_index__,colume_index__)=...
                    cov(colume_index__,rank_index__);
            end
            cov(rank_index__,rank_index__)=1+1e-6; % stabilize
            for colume_index__=rank_index__+1:x_number
                temp__=X_dis_sq(rank_index__,colume_index__,:);
                cov(rank_index__,colume_index__)=...
                    exp(-temp__(:)'*theta(:));
            end
        end
        
        % F is base funcion fval of data_point_x
%         fval_reg=ones(x_number,1); % zero
        fval_reg=[ones(x_number,1),X]; % linear
        
        % coefficient calculation
        inv_cov=inv(cov);
        beta=(fval_reg'*inv_cov*fval_reg)\fval_reg'*inv_cov*Y;
        sigma_sq=(Y-fval_reg*beta)'*inv_cov*(Y-fval_reg*beta)/x_number;
    end
    function [predict_fval,predict_variance]=interpKrigingPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            theta,beta,gama,sigma_sq,...
            inv_covariance,fval_reg,FTRF,predict_x)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        % Copyright 2022 Adel
        %
        predict_x=predict_x(:)';
        
        % normalize data
        predict_x=(predict_x-aver_X)./stdD_X;
        
        % predict value
        predict_cov__=exp(-(predict_x-X_nomlz).^2*theta(:));
        
%         predict_fval_reg__=1; % zero
        predict_fval_reg__=[1,predict_x]; % linear
        predict_fval=predict_fval_reg__*beta+predict_cov__'*gama;
        
        % predict variance
        u__=fval_reg'*inv_covariance*predict_cov__;
        predict_variance=sigma_sq*...
            (1+u__'/FTRF*u__+...
            -predict_cov__'*inv_covariance*predict_cov__);
        
        % normalize data
        predict_fval=predict_fval*stdD_Y+aver_Y;
        predict_variance=predict_variance*stdD_Y*stdD_Y;
    end
end

function interpVisualize...
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

if size(low_bou,1) ~= size(low_bou,1)
    error('interpolationRadialBasisVisualize: boundary incorrect');
end
if size(low_bou,1) > 2
    error('interpolationRadialBasisVisualize: dimension large than two');
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
    view(3);
end
end

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
