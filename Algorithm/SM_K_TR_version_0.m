clc;
clear;
close all hidden;

data_library_name='optimalSurrogateKriging_result.txt';
object_function=@func_G2_GP;
A=[];
B=[];
Aeq=[];
Beq=[];
low_bou=[-2;-2];
up_bou=[2;2];
nonlcon_function=[];

delete(data_library_name);
delete('result_total.txt');
[x_best,fval_best,NFE,output]=optimalSurrogateKriging...
    (object_function,length(low_bou),A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function);
disp(fval_best);
disp(NFE);

figure(1);
plot(output.result_fval_best);

figure(2);
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou,up_bou);
scatter3(x_list(:,1),x_list(:,2),fval_list);
xlabel('X');
ylabel('Y');
zlabel('Z');

function [x_best,fval_best,NFE,output]=optimalSurrogateKriging...
    (object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,NFE_max)
% surrogate model optimal method version 3
% compare y_predict and y_real to decide whether change range
% use kriging proxy model to decrease model_function call
% and ues ga to optimal
% iteration the model an optimal
% x_best is colume vector by ga
%
if nargin < 12
    NFE_max=30*variable_number;
    if nargin <10
        cheapcon_function=[];
        if nargin < 9
            nonlcon_function=[];
        end
    end
end
information_flag=1;
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');
data_library_name='optimalSurrogateKriging_result.txt';

% if do not input model_function, generate model_function
if nargin < 11 || isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end

done=0;NFE=0;iteration=0;
iteration_max=50;
torlance=1e-6;
sample_number=10*variable_number; % Latin hypercube sample count

low_bou_iter=low_bou;
up_bou_iter=up_bou;

bou_ratio=0.5;

result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);
iteration=iteration+1;

% check whether exist data
[x_list,~,~,~]=dataLibraryLoad...
    (data_library_name,low_bou_iter,up_bou_iter);

% generate latin hypercube sequence
[~,x_list_updata,~]=getLatinHypercube...
    (sample_number,variable_number,x_list,...
    low_bou,up_bou,cheapcon_function);

sample_number=2*variable_number; % Latin hypercube sample count
% main loop
while ~done
    
    % updata X_updata into data library
    [fval_list_updata,con_list_updata,coneq_list_updata]=dataLibraryUpdata...
        (data_library_name,model_function,x_list_updata);NFE=NFE+size(x_list_updata,1);
    
    % load data from data library
    [x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
        (data_library_name,low_bou_iter,up_bou_iter);
    
    % get kriging model and function
    [kriging_model_fval,kriging_model_con,kriging_model_coneq,output_K]=getKrigingModel...
        (x_list,fval_list,con_list,coneq_list);
    object_function=output_K.object_function;
    nonlcon_function=output_K.nonlcon_function;
    
    if information_flag
        interpolationKrigingVisualize(kriging_model_fval,low_bou_iter,up_bou_iter);
    end
    
    % updata data, x_potential
    [x_potential,fval_potential_predict]=findMinMSP...
        (object_function,variable_number,A,B,Aeq,Beq,low_bou_iter,up_bou_iter,nonlcon_function,...
        cheapcon_function);
    [fval_potential,~,~]=dataLibraryUpdata(data_library_name,model_function,x_potential');NFE=NFE+1;
    
    % find best result to record
    [x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
        (data_library_name,low_bou,up_bou);
    [x_best,fval_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        A,B,Aeq,Beq,cheapcon_function);
    
    result_x_best(iteration,:)=x_best';
    result_fval_best(iteration,:)=fval_best;
    iteration=iteration+1;
    
    if iteration > iteration_max ||...
            max(up_bou_iter-low_bou_iter) < torlance*1e4 ||...
            (iteration > 2 && ...
            abs((fval_potential-fval_potential_old)/fval_potential_old) < torlance) ||...
            NFE > NFE_max
        done=1;
    end
    
    % trust area updata guideline
    if iteration > 2
        improve=(fval_potential_old-fval_potential)/(fval_potential_old-fval_potential_predict);
        if improve <= 0.25
            bou_ratio=0.25;
        elseif improve <= 0.75
            bou_ratio=1;
        else
            bou_ratio=2;
        end
    end
    
    % updata trust area new kriging area by x_best
    interval_bou=up_bou_iter-low_bou_iter;
    interval_bou=interval_bou*bou_ratio;
    
    low_bou_iter=x_best-interval_bou/2;
    index=find(low_bou_iter < low_bou);
    low_bou_iter(index)=low_bou(index);
    
    up_bou_iter=x_best+interval_bou/2;
    index=find(up_bou_iter > up_bou);
    up_bou_iter(index)=up_bou(index);
    
    % check whether exist data
    [x_list_exist,~,~,~]=dataLibraryLoad...
        (data_library_name,low_bou_iter,up_bou_iter);
    
    % generate latin hypercube sequence ESLHS
    [~,x_list_updata,~]=getLatinHypercube...
        (sample_number,variable_number,x_list_exist,...
        low_bou_iter,up_bou_iter,cheapcon_function);
    
    x_potential_old=x_potential;
    fval_potential_old=fval_potential;
    fval_best_old=fval_best;
end

result_x_best=result_x_best(1:iteration-1);
result_fval_best=result_fval_best(1:iteration-1);

output.result_x_best=result_x_best;
output.result_fval_best=result_fval_best;

    function [fval,con,coneq]=modelFunction(x,object_function,nonlcon_function)
        % model function, concertrate favl, con, coneq into one function
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

function [x_best,fval_best]=findMinMSP...
    (object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,cheapcon_function,torlance)
% find min fval use MSP guideline
% MSP: object_funtion is object_function (generate by surrogate model)
% nonlcon_function generate by surrogate model
% use ga as optimal method
%
if ~isempty(nonlcon_function) || ~isempty(cheapcon_function)
    totalcon_function=@(x) totalconFunction(x,nonlcon_function,cheapcon_function);
else
    totalcon_function=[];
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
ga_option=optimoptions('ga','FunctionTolerance', torlance*1e6,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',50,'InitialPopulationMatrix',population_matrix);
[x_best,fval_best,~,~]=ga...
    (object_function,variable_number,A,B,Aeq,Beq,low_bou',up_bou',totalcon_function,ga_option);
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
function [x_best,fval_best]=findMinLCB(object_var_function,...
    object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,cheapcon_function,torlance)
% find min fval use LCB guideline
% LCB: object_funtion is ...
% object_function (generate by surrogate model)...
% subtract sqrt(object_var_function(x)) (generate by surrogate model)
% nonlcon_function generate by surrogate model
% use ga as optimal method
%
if ~isempty(nonlcon_function) || ~isempty(cheapcon_function)
    totalcon_function=@(x) totalconFunction(x,nonlcon_function,cheapcon_function);
else
    totalcon_function=[];
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
ga_option=optimoptions('ga','FunctionTolerance', torlance*1e6,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',50,'InitialPopulationMatrix',population_matrix,...
    'display','none');

lcb_function=@(x) object_function(x)-sqrt(object_var_function(x));
[x_best,fval_best,~,~]=ga...
    (lcb_function,variable_number,A,B,Aeq,Beq,low_bou',up_bou',totalcon_function,ga_option);

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
function [x_best,fval_best]=findMinEI(fval_list,object_var_function,nonlcon_var_function,...
    object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,cheapcon_function,torlance)
% find min fval use EI guideline
% EI: object_funtion is ei_function
% nonlcon_function generate by surrogate model
% use ga as optimal method
%

if ~isempty(nonlcon_function) || ~isempty(cheapcon_function)
    totalcon_function=@(x) totalconFunction(x,nonlcon_function,cheapcon_function);
else
    totalcon_function=[];
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
ga_option=optimoptions('ga','FunctionTolerance', torlance*1e6,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',50,'InitialPopulationMatrix',population_matrix,...
    'display','none');
fval_min=min(fval_list);
ei_function=@(x) EIFunction...
    (x,object_function,nonlcon_function,object_var_function,nonlcon_var_function,fval_min);
[x_best,fval_best,~,~]=ga...
    (ei_function,variable_number,A,B,Aeq,Beq,low_bou',up_bou',totalcon_function,ga_option);

x_best=x_best';

    function fval=EIFunction...
            (x,object_function,nonlcon_function,object_var_function,nonlcon_var_function,fval_min)
        fval=0;
        fval_predict=object_function(x);
        fval_var=object_var_function(x);
        if fval_var > 0
            normal_x=(fval_min-fval_predict)/sqrt(fval_var);
            fval=(fval_min-fval_predict)*normcdf(normal_x)+fval_var*normpdf(normal_x);
        end
        if ~isempty(nonlcon_function) && ~isempty(nonlcon_var_function)
            [con_predict_list,coneq_predict_list]=nonlcon_function(x);
            [con_var_list,coneq_var_list]=nonlcon_var_function(x);
            % process con probability
            for con_index__=1:size(con_var_list,1)
                con_predict=con_predict_list(con_index__);
                con_var=con_var_list(con_index__);
                if con_var > 0
                    fval=fval*normpdf(-con_predict/sqrt(con_var));
                else
                    if con_predict > 0
                        fval=fval*0;
                    end
                end
            end
            % process coneq probability
            for coneq_index__=1:size(coneq_var_list,1)
                coneq_predict=coneq_predict_list(coneq_index__);
                coneq_var=coneq_var_list(coneq_index__);
                if coneq_var > 0
                    fval=fval*normpdf(-coneq_predict/sqrt(coneq_var));
                else
                    if coneq_predict > 0
                        fval=fval*0;
                    end
                end
            end
        end
        fval=-fval;
    end
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
function [x_best,fval_best]=findMinRaw...
    (x_list,fval_list,con_list,coneq_list,...
    A,B,Aeq,Beq,cheapcon_function)
% find min fval in raw data
% x_list, rank is variable
% con_list, rank is con
% coneq_list, rank is coneq
% function will find min fval in con==0
% if there was not feasible x, find min consum
%
con_sum_list=zeros(size(x_list,1),1);
% add expendsive con
if ~isempty(con_list)
    con_sum_list=con_sum_list+sum(max(con_list,0),2);
end
if ~isempty(coneq_list)
    con_sum_list=con_sum_list+sum(max(coneq_list,0),2);
end

% add cheap con
for x_index=1:size(x_list,1)
    if ~isempty(A)
        if ~isempty(B)
            lincon=A*x_list(x_index,:)'-B;
        else
            lincon=A*x_list(x_index,:)';
        end
        con_sum_list(x_index)=con_sum_list(x_index)+sum(max(lincon,0));
    end
    if ~isempty(Aeq)
        if ~isempty(Beq)
            linconeq=Aeq*x_list(x_index,:)'-Beq;
        else
            linconeq=Aeq*x_list(x_index,:)';
        end
        con_sum_list(x_index)=con_sum_list(x_index)+sum(max(linconeq,0));
    end
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x_list(x_index,:)');
        con_sum_list(x_index)=con_sum_list(x_index)+...
            sum(max(con,0))+sum(max(coneq,0));
    end
end

index=find(con_sum_list==0);
if ~isempty(index)
    % feasible x
    x_list=x_list(index,:);
    fval_list=fval_list(index);
    
    % min fval
    [fval_best,index_best]=min(fval_list);
    x_best=x_list(index_best,:)';
else
    % min consum
    [fval_best,index_best]=min(con_sum_list);
    x_best=x_list(index_best,:)';
end
end

function [kriging_model_fval,kriging_model_con,kriging_model_coneq,output]=getKrigingModel...
    (x_list,fval_list,con_list,coneq_list,...
    kriging_model_fval,kriging_model_con,kriging_model_coneq)
% base on library_data to create kriging model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
% if input kriging_model_fval,kriging_model_con,kriging_model_coneq
% function will consider as updata model
% x_list,fval_list,con_list,coneq_list will regard as updata data
%
if nargin < 7
    kriging_model_coneq=[];
    if nargin < 6
        kriging_model_con=[];
        if nargin < 5
            kriging_model_fval=[];
        end
    end
end
[~,variable_number]=size(x_list);
% generate surrogate model
% base function??which can be a p x one function matrix
base_functon1=@(x) 1;
base_function_list={base_functon1};

theta_fval=ones(variable_number,1)*0.5;
if isempty(kriging_model_fval)
    kriging_model_fval=interpolationKrigingPreModel...
        (x_list,fval_list,theta_fval,base_function_list);
else
    kriging_model_fval=interpolationKrigingUpdata...
        (kriging_model_fval,x_list,fval_list);
end

if ~isempty(con_list)
    if isempty(kriging_model_fval)
        theta_con=ones(variable_number,1)*0.5;
        kriging_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
            'base_function_list',[],'F',[],'R',[],...
            'theta',[],'beta',[],'gama',[],'sigma_sq',[],...
            'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[]);
        kriging_model_con=repmat(kriging_model_con,[size(con_list,2),1]);
        for con_index=1:size(con_list,2)
            kriging_model_con(con_index)=interpolationKrigingPreModel...
                (x_list,con_list(:,con_index),theta_con,base_function_list);
        end
    else
        for con_index=1:size(con_list,2)
            kriging_model_con(con_index)=interpolationKrigingUpdata...
                (kriging_model_con(con_index),x_list,con_list(:,con_index));
        end
    end
end

if ~isempty(coneq_list)
    if isempty(kriging_model_coneq)
        theta_coneq=ones(variable_number,1)*0.5;
        kriging_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
            'base_function_list',[],'F',[],'R',[],...
            'theta',[],'beta',[],'gama',[],'sigma_sq',[],...
            'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[]);
        kriging_model_coneq=repmat(kriging_model_coneq,[size(coneq_list,2),1]);
        for coneq_index=1:size(coneq_list,2)
            kriging_model_coneq(coneq_index)=interpolationKrigingPreModel...
                (x_list,coneq_list(:,coneq_index),theta_coneq,base_function_list);
        end
    else
        for coneq_index=1:size(coneq_list,2)
            kriging_model_coneq(coneq_index)=interpolationKrigingUpdata...
                (kriging_model_coneq(coneq_index),x_new_list,coneq_list(:,coneq_index));
        end
    end
end

object_function=@(predict_x) objectFunctionSurrogate(predict_x,kriging_model_fval);
object_var_function=@(predict_x) objectVarFunctionSurrogate(predict_x,kriging_model_fval);
if isempty(kriging_model_con) && isempty(kriging_model_coneq)
    nonlcon_function=[];
    nonlcon_var_function=[];
else
    nonlcon_function=@(predict_x) nonlconFunctionSurrogate(predict_x,kriging_model_con,kriging_model_coneq);
    nonlcon_var_function=@(predict_x) nonlconVarFunctionSurrogate(predict_x,kriging_model_con,kriging_model_coneq);
end

output.object_function=object_function;
output.object_var_function=object_var_function;
output.nonlcon_function=nonlcon_function;
output.nonlcon_var_function=nonlcon_var_function;
output.x_list=x_list;
output.fval_list=fval_list;
output.con_list=con_list;
output.coneq_list=coneq_list;


    function fval=objectFunctionSurrogate...
            (predict_x,kriging_model_fval)
        [fval,~]=interpolationKrigingPredictor(kriging_model_fval,predict_x);
    end
    function fval_var=objectVarFunctionSurrogate...
            (predict_x,kriging_model_fval)
        [~,fval_var]=interpolationKrigingPredictor(kriging_model_fval,predict_x);
    end
    function [con,coneq]=nonlconFunctionSurrogate...
            (predict_x,kriging_model_con,kriging_model_coneq)
        if isempty(kriging_model_con)
            con=[];
        else
            con=zeros(length(kriging_model_con),1);
            for con_ind=1:length(kriging_model_con)
                con(con_ind)=interpolationKrigingPredictor...
                    (kriging_model_con(con_ind),predict_x);
            end
        end
        if isempty(kriging_model_coneq)
            coneq=[];
        else
            coneq=zeros(length(kriging_model_coneq),1);
            for coneq_ind=1:length(kriging_model_con)
                coneq(coneq_ind)=interpolationKrigingPredictor...
                    (kriging_model_coneq(coneq_ind),predict_x);
            end
        end
    end
    function [con_var,coneq_var]=nonlconVarFunctionSurrogate...
            (predict_x,kriging_model_con,kriging_model_coneq)
        if isempty(kriging_model_con)
            con_var=[];
        else
            con_var=zeros(length(kriging_model_con),1);
            for con_ind=1:length(kriging_model_con)
                [~,con_var(con_ind)]=interpolationKrigingPredictor...
                    (kriging_model_con(con_ind),predict_x);
            end
        end
        if isempty(kriging_model_coneq)
            coneq_var=[];
        else
            coneq_var=zeros(length(kriging_model_coneq),1);
            for coneq_ind=1:length(kriging_model_con)
                [~,coneq_var(coneq_ind)]=interpolationKrigingPredictor...
                    (kriging_model_coneq(coneq_ind),predict_x);
            end
        end
    end
end
function kriging_model=interpolationKrigingPreModel...
    (X,Y,theta,base_function_list)
% prepare model, optimal theta, load result and calculation parameter
% output is a kriging model, include X, Y, base_function_list
% theta beta gama sigma_sq is normalizede, so predict y is normalizede
% X, Y is real data
%
% Copyright 2022 Adel
%
if nargin < 4
    base_functon1=@(x) 1;
    base_function_list={base_functon1};
    if nargin < 2
        error('interpolationKrigingPreModel: lack Y');
    end
end
x_number=size(X,1);
variable_number=size(X,2);
if nargin < 3
    theta=ones(variable_number,1)*0.5;
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
X_normalize = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_normalize = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

% optimal
options = optimoptions(@fmincon,'Display','off');
lb=1e-1*ones(variable_number,1);ub=20*ones(variable_number,1);
object_function=@(theta) objectFunction(X_normalize,Y_normalize,theta,base_function_list);
theta=fmincon(object_function,theta,[],[],[],[],lb,ub,[],options);

[beta,gama,sigma_sq,F,R]=...
    interpolationKriging(X_normalize,Y_normalize,theta,base_function_list);

kriging_model.X=X;
kriging_model.Y=Y;
kriging_model.X_normalize=X_normalize;
kriging_model.Y_normalize=Y_normalize;
kriging_model.base_function_list=base_function_list;
kriging_model.F=F;
kriging_model.R=R;
kriging_model.theta=theta;
kriging_model.beta=beta;
kriging_model.gama=gama;
kriging_model.sigma_sq=sigma_sq;
kriging_model.aver_X=aver_X;
kriging_model.stdD_X=stdD_X;
kriging_model.aver_Y=aver_Y;
kriging_model.stdD_Y=stdD_Y;

    function sigma_sq=objectFunction(X,Y,theta,base_function_list)
        [~,~,sigma_sq,~,~]=interpolationKriging(X,Y,theta,base_function_list);
    end
end
function [predict_y,predict_variance]=interpolationKrigingPredictor...
    (kriging_model,predict_x)
% kriging interpolation predict function
% input predict_x and kriging model
% predict_x is row vector
% output the predict value
%
% Copyright 2022 Adel
%
if nargin < 2
    error('interpolationKrigingPredictor: lack predict_x')
end
if size(predict_x,1) > 1
    predict_x=predict_x';
end

% giving value
X_normalize=kriging_model.X_normalize;
base_function_list=kriging_model.base_function_list;
theta=kriging_model.theta;
beta=kriging_model.beta;
gama=kriging_model.gama;
sigma_sq=kriging_model.sigma_sq;
F=kriging_model.F;
R=kriging_model.R;
aver_X=kriging_model.aver_X;
stdD_X=kriging_model.stdD_X;
aver_Y=kriging_model.aver_Y;
stdD_Y=kriging_model.stdD_Y;

% normalize data
predict_x=(predict_x-aver_X)./stdD_X;

% predict value
base_function_number=size(base_function_list,1);
predict_covarience=exp(-(predict_x-X_normalize).^2*theta);
predict_F=zeros(base_function_number,1);
for base_function_index=1:base_function_number
    base_function=base_function_list{base_function_index};
    predict_F(base_function_index,1)=base_function(predict_x);
end
predict_F=predict_F';
predict_y=predict_F'*beta+predict_covarience'*gama;

% predict variance
u=F'/R*predict_covarience-predict_F;
predict_variance=sigma_sq*(1+u'/(F'/R*F)*u-predict_covarience'/R*predict_covarience);

% normalize data
predict_y=predict_y*stdD_Y+aver_Y;
predict_variance=predict_variance*stdD_Y*stdD_Y;
end
function kriging_model=interpolationKrigingUpdata...
    (kriging_model,X_new,Y_new)
% add new point in existing kriging model
% x_new is row vector, can be more than one
% y_new is row vector, can be more than one
%

% giving value
X=kriging_model.X;
Y=kriging_model.Y;
base_function_list=kriging_model.base_function_list;
theta=kriging_model.theta;

% add new data
X=[X;X_new];
Y=[Y;Y_new];

x_number=size(X,1);

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__ = find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__) = 1; end
index__ = find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__) = 1; end
X_normalize = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_normalize = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

% updata model
[beta,gama,sigma_sq,F,R]=...
    interpolationKriging(X_normalize,Y_normalize,theta,base_function_list);

kriging_model.X=X;
kriging_model.Y=Y;
kriging_model.X_normalize=X_normalize;
kriging_model.Y_normalize=Y_normalize;
kriging_model.beta=beta;
kriging_model.gama=gama;
kriging_model.sigma_sq=sigma_sq;
kriging_model.F=F;
kriging_model.R=R;
kriging_model.aver_X=aver_X;
kriging_model.stdD_X=stdD_X;
kriging_model.aver_Y=aver_Y;
kriging_model.stdD_Y=stdD_Y;

end
function interpolationKrigingVisualize...
    (kriging_model,low_bou,up_bou,figure_handle)
% visualization polynamial respond surface model
% figrue is 100
%
if nargin < 4
    figure_handle=figure(100);
end
if size(low_bou,1) ~= size(low_bou,1)
    error('interpolationKrigingVisualize: boundary incorrect');
end
if size(low_bou,1) > 2
    error('interpolationKrigingVisualize: dimension large than two');
end
axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end

x_list=kriging_model.X;
y_list=kriging_model.Y;

grid_number=100;
d_bou=(up_bou-low_bou)/grid_number;
predict_result=zeros(grid_number+1);
[X,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
for x_index=1:grid_number+1
    for y_index=1:grid_number+1
        predict_x=([x_index,y_index]-1).*d_bou'+low_bou';
        [predict_result(y_index,x_index),~]=interpolationKrigingPredictor...
            (kriging_model,predict_x);
    end
end
surf(axes_handle,X,Y,predict_result,'FaceAlpha',0.5,'EdgeColor','none');
line(axes_handle,x_list(:,1),x_list(:,2),y_list,'Marker','o','LineStyle','none');
xlabel('X');
ylabel('Y');
zlabel('Z');
end
function [beta,gama,sigma_sq,F,R]=interpolationKriging...
    (X,Y,theta,base_function_list)
% total riging interpolation function
% input X, Y as initial data, theta and base function
% output beta and gama, xita
%
% Copyright 2022 Adel
%
x_number=size(X,1);

% covariance matrix
R=zeros(x_number);
for rank=1:x_number
    for colume=1:rank-1
        R(rank,colume)=R(colume,rank);
    end
    for colume=rank:x_number
        R(rank,colume)=exp(-(X(rank,:)-X(colume,:)).^2*theta);
    end
end
base_function_number=size(base_function_list,1);

% F is base funcion fval of data_point_x
F=zeros(base_function_number,x_number);
for base_function_index=1:base_function_number
    base_function=base_function_list{base_function_index};
    for x_index=1:x_number
        F(base_function_index,x_index)=base_function(X(x_index,:));
    end
end
F=F';

if rcond(R)<1e-16
    disp('rcond error');
end
% coefficient calculation
beta=(F'/R*F)\F'/R*Y;
gama=R\(Y-F*beta);
sigma_sq=(Y-F*beta)'/R*(Y-F*beta)/x_number;
end

function [fval_list,con_list,coneq_list]=dataLibraryUpdata...
    (data_library_name,model_function,x_list)
% updata data library
%
fval_list=[];
con_list=[];
coneq_list=[];
[x_number,variable_number]=size(x_list);

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
                con_list=[con_list;...
                    data(5+variable_number:4+variable_number+con_number)];
                coneq_list=[coneq_list;...
                    data(6+variable_number+con_number:5+variable_number+con_number+coneq_number)];
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

function [X,X_new,distance_min_normalize]=getLatinHypercube...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou,cheapcon_function)
% generate sample sequence latin hypercube
% iteration optimal method is used
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
            if nargin < 2
                error('getLatinHypercube: lack variable_number');
            end
        end
        low_bou=zeros(variable_number,1);
        up_bou=ones(variable_number,1);
    end
end
iteration_max=20*variable_number;

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
    X_exist_norml=(X_exist-low_bou')./(up_bou'-low_bou');
else
    X_exist_norml=[];
end

% check x_new_number
x_new_number=sample_number-size(X_exist,1);
if x_new_number < 0
    X=X_exist;
    X_new=[];
    distance_min_normalize=getMinDistance(X_exist_norml);
    return;
end

low_bou_norml=zeros(variable_number,1);
up_bou_norml=ones(variable_number,1);

% initialize X_new by sample with constraint function
[X_new,distance_min_normalize]=getLatinHypercubeInitial...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou,cheapcon_function);

% x is nomalize, so constraint function should change
if ~isempty(cheapcon_function)
    cheapcon_function=@(x) max(sample_number*cheapcon_function(x.*(up_bou-low_bou)+low_bou)+1,0);
end

iteration=0;
while iteration < iteration_max
    %     scatter(X_new(:,1),X_new(:,2));
    %     bou=[low_bou_norml,up_bou_norml]';
    %     axis(bou(:));axis equal
    %     radius=1/6;
    %     hold on;
    %     rectangle('Position',[-radius+0.5,-radius+0.5,2*radius,2*radius],'Curvature',[1 1])
    %     hold off;
    
    % change each x place by newton methods
    fval_list=zeros(x_new_number,1);
    gradient_list=zeros(variable_number,x_new_number);
    for X_new_index=1:x_new_number
        % get gradient
        [fval_list(X_new_index,1),gradient_list(:,X_new_index),~]=objectFunctionXPlace...
            (X_new(X_new_index,:)',[X_new(1:X_new_index-1,:);X_new(X_new_index+1:end,:);X_exist_norml],...
            variable_number,low_bou_norml,up_bou_norml,cheapcon_function);
    end
    
    fval_list=fval_list/max(fval_list);
    for X_new_index=1:x_new_number
        x=X_new(X_new_index,:)'-...
            fval_list(X_new_index,1)*gradient_list(:,X_new_index)/norm(gradient_list(:,X_new_index),2)*...
            distance_min_normalize*(1-iteration/iteration_max);
        for variable_index=1:variable_number
            if x(variable_index,1) > 1
                x(variable_index,1)=1;
            end
            if x(variable_index,1) < 0
                x(variable_index,1)=0;
            end
        end
        X_new(X_new_index,:)=x';
    end
    
    iteration=iteration+1;
end
distance_min_normalize=getMinDistance([X_new;X_exist_norml]);
X_new=X_new.*(up_bou'-low_bou')+low_bou';
X=[X_new;X_exist];

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
        dimension=size(low_bou,1);
        
        switch dimension
            case 1
                d_bou=(up_bou-low_bou)/grid_number;
                X__=low_bou:d_bou:up_bou;
                fval__=zeros(grid_number+1,1);
                for x_index=1:(grid_number+1)
                    fval__(x_index)=object_function(X__(x_index));
                end
                line(axes_handle,X__,fval__);
                xlabel('X');
                ylabel('value');
                
            case 2
                d_bou=(up_bou-low_bou)/grid_number;
                [X__,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
                fval__=zeros(grid_number+1);
                for x_index=1:grid_number+1
                    for y_index=1:grid_number+1
                        predict_x=([x_index;y_index]-1).*d_bou+low_bou;
                        fval__(x_index,y_index)=object_function(predict_x);
                    end
                end
                fval__(find(fval__ > Y_max))=Y_max;
                fval__(find(fval__ < Y_min))=Y_min;
                surf(axes_handle,X__',Y',fval__,'FaceAlpha',0.5,'EdgeColor','none');
                xlabel('X');
                ylabel('Y');
                zlabel('value');
        end
        
    end
    function [fval,gradient,hessian]=objectFunctionXPlace...
            (x,X_surplus,variable_number,low_bou,up_bou,cheapcon_function)
        % function describe distance between X and X_supply
        % X__ is colume vector and X_supply__ is matrix which is num-1 x var
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
            
            fval__(:,2)=differ_function(x);
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
    function [X_new_norml,distance_min_normalize__]=getLatinHypercubeInitial...
            (sample_number,variable_number,X_exist,...
            low_bou,up_bou,cheapcon_function)
        % generate sample sequence latin hypercube
        % election sequential method is used
        % sample number is total point in area
        % default low_bou is 0, up_bou is 1, cheapcon_function is []
        % low_bou and up_bou is colume vector
        % x in x_exist_list, x_list, supply_x_list is row vector
        % x_exist_list should meet bou
        %
        % Copyright 2022 Adel
        %
        
        iteration_max__=10*variable_number;
        x_new_number__=sample_number-size(X_exist,1);
        
        % get quasi-feasible point
        iteration__=0;
        x_initial_number=100*x_new_number__;
        X_supply_initial_normalize=rand(x_initial_number,variable_number);
        if ~isempty(cheapcon_function)
            x_index=1;
            while x_index <= size(X_supply_initial_normalize,1)
                x_supply=X_supply_initial_normalize(x_index,:).*(up_bou'-low_bou')+low_bou';
                if cheapcon_function(x_supply) > 0
                    X_supply_initial_normalize(x_index,:)=[];
                else
                    x_index=x_index+1;
                end
            end
            % check if have enough X_supply_normalize
            while size(X_supply_initial_normalize,1) < 10*x_new_number__ &&...
                    iteration__ < iteration_max__
                
                X_supply_initial_normalize_new=rand(x_initial_number,variable_number);
                x_index=1;
                while x_index <= size(X_supply_initial_normalize_new,1)
                    x_supply=X_supply_initial_normalize_new(x_index,:).*(up_bou'-low_bou')+low_bou';
                    if cheapcon_function(x_supply) > 0
                        X_supply_initial_normalize_new(x_index,:)=[];
                    else
                        x_index=x_index+1;
                    end
                end
                
                % add more point
                X_supply_initial_normalize=[X_supply_initial_normalize;...
                    X_supply_initial_normalize_new];
                iteration__=iteration__+1;
            end
        end
        X_supply_quasi_normalize=X_supply_initial_normalize;
        
        % iterate and get final x_supply_list
        iteration__=0;
        x_supply_quasi_number=size(X_supply_quasi_normalize,1);
        if ~isempty(X_exist)
            x_exist_list_normalize=(X_exist-low_bou')./(up_bou'-low_bou');
        else
            x_exist_list_normalize=[];
        end
        distance_min_normalize__=0;
        X_new_norml=[];
        while iteration__ <= iteration_max__
            % random select x_new_number X to X_trial_normalize
            x_select_index=randperm(x_supply_quasi_number,x_new_number__);
            X_trial_normalize=[X_supply_quasi_normalize(x_select_index,:);x_exist_list_normalize];
            
            % get distance min itertion X_
            distance_min_iteration=getMinDistance(X_trial_normalize);
            
            % if distance_min_iteration is large than last time
            if distance_min_iteration > distance_min_normalize__
                distance_min_normalize__=distance_min_iteration;
                X_new_norml=X_supply_quasi_normalize(x_select_index,:);
            end
            
            iteration__=iteration__+1;
        end
    end
end