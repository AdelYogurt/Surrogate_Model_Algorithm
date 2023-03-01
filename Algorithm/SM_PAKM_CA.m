
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

% variable_number=2;
% object_function=@(x) benchmark.singleWeiObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[0,0];
% up_bou=[3.7,4];
% nonlcon_function=@(x) benchmark.singleWeiNonlcon(x);
% cheapcon_function=[];

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
[x_best,fval_best,NFE,output]=optimalSurrogatePAKMCA...
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
function [x_best,fval_best,NFE,output]=optimalSurrogatePAKMCA...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance)
% Parallel Adaptive Kriging Method with Constraint Aggregation
%
% x_list is x_number x variable_number matrix
% both nonlcon_function and cheapcon_function format is [con,coneq]
% model_function should output fval, format is [fval,con,coneq]
% con or coneq can be colume vector if there was more than one constrain
%
% referance: LONG T, WEI Z, SHI R, et al. Parallel Adaptive Kriging Method
% with Constraint Aggregation for Expensive Black-Box Optimization Problems
% [J]. AIAA Journal, 2021, 59(9): 3465-79.
%
% Copyright Adel 2023.2
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
CONVERGENCE_JUDGMENT_FLAG=0; % whether judgment convergence

if isempty(iteration_max)
    iteration_max=100;
end

% hyper parameter
if variable_number < 10
    sample_number_initial=min((variable_number+1)*(variable_number+2)/2,5*variable_number);
else
    sample_number_initial=variable_number+1;
end
if variable_number <= 5
    sample_number_iteration=2;
else
    sample_number_iteration=3;
end
rou=4;
rou_min=1;
rou_max=64;
rou_decrease=0.5;
rou_increase=2;

% max fval when normalize fval, con, coneq
nomlz_fval=1;

protect_range=1e-4;

data_library_name='optimal_data_library';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;

result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);

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

kriging_model_fval=[];
kriging_model_con=[];
kriging_model_coneq=[];
kriging_model_KS=[];
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
    fval_nomlz_list=fval_list./fval_max*nomlz_fval;
    if ~isempty(con_list)
        con_max_list=max(abs(con_list),[],1);
        con_nomlz_list=con_list./con_max_list*nomlz_fval;
    else
        con_max_list=[];
        con_nomlz_list=[];
    end
    if ~isempty(coneq_list)
        coneq_max_list=max(abs(coneq_list),[],1);
        coneq_nomlz_list=coneq_list./coneq_max_list*nomlz_fval;
    else
        coneq_max_list=[];
        coneq_nomlz_list=[];
    end

    % step 4
    [kriging_model_fval,kriging_model_con,kriging_model_coneq,output_kriging]=getKrigingModel...
        (x_list,fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
        kriging_model_fval,kriging_model_con,kriging_model_coneq);
    object_function_surrogate=output_kriging.object_function_surrogate;
    nonlcon_function_surrogate=output_kriging.nonlcon_function_surrogate;

    % step 5
    % MSP guideline to obtain x_adapt
%     if expensive_nonlcon_flag
%         nonlcon_torlance_surrogate=min(max(con_nomlz_list,[],2))*(1-NFE/NFE_max)^2;
%     else
%         nonlcon_torlance_surrogate=0;
%     end
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);
    [x_best_pred,~,exitflag,~]=findMinMSP...
        (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function_surrogate,...
        cheapcon_function,nonlcon_torlance,x_best);

    if exitflag == -2
        % optimal feasiblilty if do not exist feasible point
        object_nonlcon_function_surrogate=@(x) objectNonlconFunctionSurrogate(x,nonlcon_function_surrogate);
        [x_best_pred,~,exitflag,~]=findMinMSP...
            (object_nonlcon_function_surrogate,variable_number,low_bou,up_bou,[],...
            cheapcon_function,nonlcon_torlance,x_best);
    end

    % check x_potential if exist in data library
    % if not, updata data libraray
    [x_best_pred,fval_best_pred,con_best_pred,coneq_best_pred,NFE_p,repeat_index]=dataLibraryUpdataProtect...
        (data_library_name,model_function,x_best_pred,x_list,low_bou,up_bou,protect_range);NFE=NFE+NFE_p;
    x_list=[x_list;x_best_pred];
    fval_list=[fval_list;fval_best_pred];
    if ~isempty(con_list)
        con_list=[con_list;con_best_pred];
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_best_pred];
    end

    % normalization fval updata
    if ~isempty(fval_best_pred)
        fval_potential_nomlz=fval_best_pred/fval_max*nomlz_fval;
        fval_nomlz_list=[fval_nomlz_list;fval_potential_nomlz];
    end
    if ~isempty(con_best_pred)
        con_potential_nomlz=(con_best_pred./con_max_list)*nomlz_fval;
        con_nomlz_list=[con_nomlz_list;con_potential_nomlz];
    end
    if ~isempty(coneq_best_pred)
        coneq_potential_nomlz=(coneq_best_pred./coneq_max_list)*nomlz_fval;
        coneq_nomlz_list=[coneq_nomlz_list;coneq_potential_nomlz];
    end

    % when x_potential is exist in data library, x_potential_add will be
    % empty, this times we will use origin point data
    if isempty(x_best_pred)
        x_best_pred=x_list(repeat_index,:);
        fval_best_pred=fval_list(repeat_index,:);
        if ~isempty(con_list)
            con_best_pred=con_list(repeat_index,:);
        else
            con_best_pred=[];
        end
        if ~isempty(coneq_list)
            coneq_best_pred=coneq_list(repeat_index,:);
        else
            coneq_best_pred=[];
        end
    end

    if DRAW_FIGURE_FLAG && variable_number < 3
        interpVisualize(kriging_model_fval,low_bou,up_bou);
        line(x_best_pred(1),x_best_pred(2),fval_best_pred/fval_max*nomlz_fval,'Marker','o','color','r','LineStyle','none')
    end

    % step 6
    % find best result to record
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function,nonlcon_torlance);

    if INFORMATION_FLAG
        fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
        fprintf('current x:          %s\n',num2str(x_best_pred));
        fprintf('current value:      %f\n',fval_best_pred);
        fprintf('current violation:  %s  %s\n',num2str(con_best_pred),num2str(coneq_best_pred));
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
                abs((fval_best_pred-fval_potential_old)/fval_potential_old) < torlance)
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

    % PCFEI Function-Based Infill Sampling Mechanism
    if ~done
        % step 1
        % construct kriging model of KS function
        % updata rou
        if (max(con_best_pred) < nonlcon_torlance)
            rou=rou*rou_increase;
        else
            rou=rou*rou_decrease;
        end
        rou=max(rou,rou_min);
        rou=min(rou,rou_max);

        % modify
%         KS_nomlz=max(con_nomlz_list,[],2);
        KS_nomlz=log(sum(exp(con_nomlz_list*rou),2))/rou;
        [kriging_model_KS,~,~,output]=getKrigingModel...
            (x_list,KS_nomlz,[],[],...
            kriging_model_KS);
        object_function_surrogate_KS=output.object_function_surrogate;

        % step 2
        % contruct EI, PF function
        object_function_EI=@(X) EIFunction(object_function_surrogate,X,fval_best/fval_max);
        object_function_PF=@(X) PFFunction(object_function_surrogate_KS,X);
        object_function_IF=@(X) IFFunction(x_best_pred,X,exp(kriging_model_fval.hyp),variable_number);

        % step 3
        % multi objective optimization to get pareto front
        object_function_PCFEI=@(x) [-object_function_EI(x),-object_function_PF(x)];
        gamultiobj_option=optimoptions('gamultiobj','Display','none');
        [x_pareto_list,fval_pareto_list,exitflag,output_gamultiobj]=gamultiobj...'
            (object_function_PCFEI,variable_number,[],[],[],[],low_bou,up_bou,[],gamultiobj_option);

        % step 4
        % base on PCFEI value to get first sample_number_iteration point
        EI_list=-fval_pareto_list(:,1);
        EI_list=EI_list/max(EI_list);
        PF_list=-fval_pareto_list(:,2);
        PF_list=PF_list/max(PF_list);
        IF_list=object_function_IF(x_pareto_list);
        IF_list=IF_list/max(IF_list);
        PCFEI_list=EI_list.*PF_list.*IF_list;
        [~,index_list]=sort(PCFEI_list);
        x_updata_list=x_pareto_list(index_list((end+1-sample_number_iteration):end),:);
    end

    x_potential_old=x_best_pred;
    fval_potential_old=fval_best_pred;
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

    function fval=EIFunction(object_function_surrogate,X,fval_min)
        % EI function
        [Fval_pred,Fval_var]=object_function_surrogate(X);
        normal_fval=(fval_min-Fval_pred)./sqrt(Fval_var);
        EI_l=(fval_min-Fval_pred).*normcdf(normal_fval);
        EI_g=Fval_var.*normpdf(normal_fval);
        fval=EI_l+EI_g;
    end

    function fval=PFFunction(object_function_surrogate,X)
        % PF function
        [Con_pred,Con_var]=object_function_surrogate(X);
        fval=normcdf(-Con_pred./sqrt(Con_var));
    end

    function fval=IFFunction(x_best,X,theta,variable_number)
        fval=zeros(size(X,1),1);
        for variable_index=1:variable_number
            fval=fval+(X(:,variable_index)-x_best(:,variable_index)').^2*theta(variable_index);
        end
        fval=1-exp(-fval);
    end
end

%% auxiliary function
function [x_best,fval_best,exitflag,output]=findMinMSP...
    (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function_surrogate,...
    cheapcon_function,nonlcon_torlance,x_best)
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
population_matrix=zeros(min(10,2*variable_number),variable_number);
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
    'PopulationSize',min(10,2*variable_number),...
    'MaxGenerations',50,'InitialPopulationMatrix',population_matrix,...
    'display','none');
[x_best,fval_best,exitflag,output]=ga...
    (object_function_surrogate,variable_number,[],[],[],[],low_bou,up_bou,constraint_function,ga_option);
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
        [con,coneq]=cheapcon_function(x_list(x_index,:));
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
    x_best=x_list(index_best,:);
    if ~isempty(con_list)
        con_best=con_list(index_best,:);
    end
    if ~isempty(coneq_list)
        coneq_best=coneq_list(index_best,:);
    end
else
    % min consum
    [~,index_best]=min(max_nonlcon_list);
    fval_best=fval_list(index_best);
    x_best=x_list(index_best,:);
    if ~isempty(con_list)
        con_best=con_list(index_best,:);
    end
    if ~isempty(coneq_list)
        coneq_best=coneq_list(index_best,:);
    end
end
end

%% surrogate model
function [kriging_model_fval,kriging_model_con,kriging_model_coneq,output]=getKrigingModel...
    (x_list,fval_list,con_list,coneq_list,...
    kriging_model_fval,kriging_model_con,kriging_model_coneq)
% base on library_data to create kriging model and function
% if input model, function will updata model
% object_function is multi fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
if size(x_list,1) ~= size(fval_list,1)
    error('getKrigingModel: x_list size no equal fval_list size')
end

if isempty(kriging_model_fval)
    [predict_function_fval,kriging_model_fval]=interpKrigingPreModel...
        (x_list,fval_list);
else
    [predict_function_fval,kriging_model_fval]=interpKrigingPreModel...
        (x_list,fval_list,kriging_model_fval.hyp);
end

if ~isempty(con_list)
    predict_function_con=cell(size(con_list,2),1);
    if size(x_list,1) ~= size(con_list,1)
        error('getKrigingModel: x_list size no equal con_list size')
    end
    if isempty(kriging_model_con)
        kriging_model_con=struct('X',[],'Y',[],...
            'fval_regression',[],'covariance',[],'inv_covariance',[],...
            'hyp',[],'beta',[],'gama',[],'sigma_sq',[],...
            'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
            'predict_function',[]);
        kriging_model_con=repmat(kriging_model_con,1,[size(con_list,2)]);
        for con_index=1:size(con_list,2)
            [predict_function_con{con_index},kriging_model_con(con_index)]=interpKrigingPreModel...
                (x_list,con_list(:,con_index));
        end
    else
        for con_index=1:size(con_list,2)
            [predict_function_con{con_index},kriging_model_con(con_index)]=interpKrigingPreModel...
                (x_list,con_list(:,con_index),kriging_model_con(con_index).hyp);
        end
    end
else
    predict_function_con=[];
    kriging_model_con=[];
end

if ~isempty(coneq_list)
    predict_function_coneq=cell(size(coneq_list,2),1);
    if size(x_list,1) ~= size(coneq_list,1)
        error('getKrigingModel: x_list size no equal coneq_list size')
    end
    if isempty(kriging_model_coneq)
        kriging_model_coneq=struct('X',[],'Y',[],...
            'fval_regression',[],'covariance',[],'inv_covariance',[],...
            'hyp',[],'beta',[],'gama',[],'sigma_sq',[],...
            'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[],...
            'predict_function',[]);
        kriging_model_coneq=repmat(kriging_model_coneq,1,[size(coneq_list,2)]);
        for coneq_index=1:size(coneq_list,2)
            [predict_function_coneq{coneq_index},kriging_model_coneq(coneq_index)]=interpKrigingPreModel...
                (x_list,coneq_list(:,coneq_index));
        end
    else
        for coneq_index=1:size(coneq_list,2)
            [predict_function_coneq{coneq_index},kriging_model_coneq(coneq_index)]=interpKrigingPreModel...
                (x_list,coneq_list(:,coneq_index),kriging_model_coneq(coneq_index).hyp);
        end
    end
else
    predict_function_coneq=[];
    kriging_model_coneq=[];
end

object_function_surrogate=@(X_predict) objectFunctionSurrogate(X_predict,predict_function_fval);
if isempty(predict_function_con) && isempty(predict_function_coneq)
    nonlcon_function_surrogate=[];
else
    nonlcon_function_surrogate=@(X_predict) nonlconFunctionSurrogate(X_predict,predict_function_con,predict_function_coneq);
end

output.object_function_surrogate=object_function_surrogate;
output.nonlcon_function_surrogate=nonlcon_function_surrogate;

    function [fval,fval_var]=objectFunctionSurrogate...
            (X_predict,predict_function_fval)
        % connect all predict favl
        %
        [fval,fval_var]=predict_function_fval(X_predict);
    end

    function [con,con_var,coneq,coneq_var]=nonlconFunctionSurrogate...
            (X_predict,predict_function_con,predict_function_coneq)
        % connect all predict con and coneq
        %
        if isempty(predict_function_con)
            con=[];
            con_var=[];
        else
            con=zeros(size(X_predict,1),length(predict_function_con));
            con_var=zeros(size(X_predict,1),length(predict_function_con));
            for con_index__=1:length(predict_function_con)
                [con(:,con_index__),con_var(:,con_index__)]=....
                    predict_function_con{con_index__}(X_predict);
            end
        end
        if isempty(predict_function_coneq)
            coneq=[];
            coneq_var=[];
        else
            coneq=zeros(size(X_predict,1),length(predict_function_coneq));
            coneq_var=zeros(size(X_predict,1),length(predict_function_coneq));
            for coneq_index__=1:length(predict_function_coneq)
                [coneq(:,coneq_index__),coneq_var(:,coneq_index__)]=...
                    predict_function_coneq{coneq_index__}(X_predict);
            end
        end
    end
end

function [predict_function,kriging_model]=interpKrigingPreModel...
    (X,Y,hyp)
% nomalization method is grassian
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
% reg_function=@(X) regZero(X);
reg_function=@(X) regLinear(X);

% calculate reg
fval_reg_nomlz=(reg_function(X)-aver_Y)./stdD_Y;

% optimal to get hyperparameter
fmincon_option=optimoptions('fmincon','Display','none',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',false);
low_bou_hyp=-3*ones(1,variable_number);
up_bou_hyp=3*ones(1,variable_number);
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
        fval_reg_pred_nomlz=(fval_reg_pred-aver_Y)./stdD_Y;
        
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
