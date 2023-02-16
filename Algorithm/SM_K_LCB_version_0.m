clc;
clear;
close all hidden;

benchmark_function=BenchmarkFunction();

data_library_name='optimal_SM_K_LCB_result.txt';

% variable_number=2;
% object_function=@(x) benchmark_function.singlePKObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-3,-3];
% up_bou=[3,3];
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

object_function=@(x) benchmark_function.singleHNObject(x);
variable_number=6;
A=[];
B=[];
Aeq=[];
Beq=[];
low_bou=zeros(1,variable_number);
up_bou=ones(1,variable_number);
nonlcon_function=[];
cheapcon_function=[];
model_function=[];

[x_best,fval_best,NFE,output]=optimalSurrogateKLCB...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function)
result_x_best=output.result_x_best;
result_fval_best=output.result_fval_best;

figure(1);
plot(result_fval_best);

function [x_best,fval_best,NFE,output]=optimalSurrogateKLCB...
    (object_function_surrogate,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,NFE_max)
% surrogate base optimal method version 0
% use LCB function to decide which point was select to updata and use ga to
% optimal.
% unsupport nonlcon.
% model_function should output fval, format is [fval,con,coneq]
%
% Copyright 2022 Adel
%
if nargin < 8
    if nargin < 7
        cheapcon_function=[];
        if nargin < 6
            nonlcon_function=[];
        end
    end
    NFE_max=20*variable_number;
end
DRAW_FIGURE_FLAG=0; % whether draw data
CONVERGENCE_JUDGMENT_FLAG=1; % whether judgment convergence

iteration_max=100;
torlance=1e-3;
sample_number_initial=5*variable_number;

w_list=[1];

data_library_name='optimal_SM_K_LCB_result.txt';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;
result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);

% if do not input model_function, generate model_function
if nargin < 11 || isempty(model_function)
    model_function=@(x) modelFunction(x,object_function_surrogate,nonlcon_function);
end

iteration=iteration+1;

% generate latin hypercube sequence
[~,x_list_updata,~]=getLatinHypercube...
    (sample_number_initial,variable_number,[],...
    low_bou,up_bou,cheapcon_function);

x_list=[];
fval_list=[];
kriging_model_fval=[];
while ~done
    % updata X_updata into data library
    [fval_list_updata,~,~]=dataLibraryUpdata...
        (data_library_name,model_function,x_list_updata);NFE=NFE+size(x_list_updata,1);

    x_list=[x_list;x_list_updata];
    fval_list=[fval_list;fval_list_updata];

    % nomalization con
    fval_max=max(abs(fval_list));
    fval_nomlz_list=fval_list/fval_max*10;
    
        % step 4
    % generate kriging model use normalization fval
    % if x_list more than 11*D-1+25, select 11*D-1+25 to construct model
    if size(x_list,1) > (11*variable_number-1+25)
        % select by nearest to x_best
        distance=sum(abs(x_best-x_list),2);
        [~,nearest_index_list]=sort(distance);
        nearest_index_list=nearest_index_list(5*variable_number);

        % rand select from remainder
        index_list=1:size(x_list,1);
        index_list(nearest_index_list)=[];
        index_list=index_list(randi(length(index_list),...
            [1,(11*variable_number-1+25-length(nearest_index_list))]));

        x_list_model=[x_list(nearest_index_list,:);x_list(index_list,:)];
        fval_nomlz_model_list=[fval_nomlz_list(nearest_index_list,:);fval_nomlz_list(index_list,:)];
    else
        x_list_model=x_list;
        fval_nomlz_model_list=fval_nomlz_list;
    end
    [kriging_model_fval,output_K]=getKrigingModel...
        (x_list_model,fval_nomlz_model_list,kriging_model_fval);
    object_function_surrogate=output_K.object_function_surrogate;
    object_function_variance=output_K.object_function_variance;
    
    % find best result to record
    [x_best,fval_best]=findMinRaw...
        (x_list,fval_list,[],[],...
        cheapcon_function);
    
    result_x_best(iteration,:)=x_best;
    result_fval_best(iteration,:)=fval_best;
    iteration=iteration+1;
    
    % LCB guideline to obtain x_adapt
    w=w_list(mod(iteration-1,length(w_list))+1);
    [x_adapt,~,LCB_function]=findMinLCB...
        (object_function_surrogate,variable_number,low_bou,up_bou,cheapcon_function,...
        object_function_variance,w);
    
    if DRAW_FIGURE_FLAG
        interpolationVisualize(kriging_model_fval,low_bou,up_bou);
        drawFunction(lcb_function,low_bou,up_bou)
    end
    
    % forced interrupt
    if iteration > iteration_max || NFE > NFE_max
        done=1;
    end
    
    % convergence judgment
    if CONVERGENCE_JUDGMENT_FLAG
        distance=sum((x_list-x_adapt).^2./(low_bou-up_bou).^2,2);
        % if the gain of change proxy is too less, stop proxy iteration
        if  min(distance) < torlance*1e-4
            done=1;
        end
    end
    
    x_list_updata=x_adapt;
end
result_x_best=result_x_best(1:iteration-1,:);
result_fval_best=result_fval_best(1:iteration-1);

output.result_x_best=result_x_best;
output.result_fval_best=result_fval_best;

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

function [x_best,fval_best,LCB_function]=findMinLCB...
    (object_function_surrogate,variable_number,low_bou,up_bou,cheapcon_function,...
    object_function_variance,w)
% find min fval use LCB guideline
% LCB: object_funtion is ...
% object_function (generate by surrogate model)...
% subtract sqrt(object_var_function(x)) (generate by surrogate model)
% nonlcon_function generate by surrogate model
% use ga as optimal method
%
if nargin < 12
    w=1;
end

% generate initial population for ga
population_matrix=zeros(2*variable_number,variable_number);
for population_index=1:2*variable_number
    x=rand(1,variable_number).*(up_bou-low_bou)+low_bou;
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x);
        while sum(~(con < 0))
            x=rand(1,variable_number).*(up_bou-low_bou)+low_bou;
            [con,coneq]=cheapcon_function(x);
        end
    end
    population_matrix(population_index,:)=x;
end

% optiaml
ga_option=optimoptions('ga','FunctionTolerance', 1e-2,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',50,'InitialPopulationMatrix',population_matrix,...
    'display','none');

LCB_function=@(x) object_function_surrogate(x)-w*sqrt(object_function_variance(x));
[x_best,fval_best,~,~]=ga...
    (LCB_function,variable_number,[],[],[],[],low_bou,up_bou,cheapcon_function,ga_option);
end

function [x_best,fval_best]=findMinRaw...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function)
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
    x_best=x_list(index_best,:);
else
    % min consum
    [fval_best,index_best]=min(con_sum_list);
    x_best=x_list(index_best,:);
end
end

function [kriging_model_fval,output]=getKrigingModel...
    (x_list,fval_list,kriging_model_fval)
% base on library_data to create kriging model and function
%
% object_function is single fval output
%
% if input kriging_model_fval,kriging_model_con,kriging_model_coneq
% function will consider as updata model
%
if size(x_list,1) ~= size(fval_list,1)
    error('getKrigingModel: x_list size no equal fval_list size')
end
% generate surrogate model
if isempty(kriging_model_fval)
    [predict_function,kriging_model_fval]=interpKrigingPreModel...
        (x_list,fval_list);
else
    [predict_function,kriging_model_fval]=interpKrigingPreModel...
        (x_list,fval_list,kriging_model_fval.theta);
end

object_function_surrogate=@(predict_x) objectFunctionSurrogate(predict_x,predict_function);
object_function_variance=@(predict_x) objectFunctionVariance(predict_x,predict_function);

output.object_function_surrogate=object_function_surrogate;
output.object_function_variance=object_function_variance;
output.x_list=x_list;
output.fval_list=fval_list;

    function fval=objectFunctionSurrogate...
            (predict_x,predict_function)
        [fval,~]=predict_function(predict_x);
    end
    function fval_var=objectFunctionVariance...
            (predict_x,predict_function)
        [~,fval_var]=predict_function(predict_x);
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
fmincon_option=optimoptions(@fmincon,'Display','off',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10);
low_bou_kriging=1e-1*ones(1,variable_number);
up_bou_kriging=20*ones(1,variable_number);
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
if size(low_bou,1) ~= size(low_bou,1)
    error('interpolationRadialBasisVisualize: boundary incorrect');
end
if size(low_bou,1) > 2
    error('interpolationRadialBasisVisualize: dimension large than two');
end

axes_handle=axes(figure_handle);

x_list=model.X;
y_list=model.Y;
predict_function=model.predict_function;

% get boundary
if isempty(low_bou)
    low_bou=min(x_list,[],1);
else
    low_bou=low_bou(:)';
end
if isempty(up_bou)
    up_bou=max(x_list,[],1);
else
    up_bou=up_bou(:)';
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
            predict_x=([x_index,y_index]-1).*d_bou+low_bou;
            [predict_result(y_index,x_index),~]=predict_function(predict_x);
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
% updata format:
% variable_number, fval_number, con_number, coneq_number
% x, fval, con, coneq
%
fval_list=[];
con_list=[];
coneq_list=[];
[x_number,variable_number]=size(x_list);

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