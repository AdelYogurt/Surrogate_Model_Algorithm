clc;
clear;
close all hidden;

data_library_name='optimalSurrogateSADEKTS_result';

variable_number=4;
object_function_H=@(x) func_H2_ROS_H(x);
object_function_L=@(x) func_H2_ROS_L(x);
A=[];
B=[];
Aeq=[];
Beq=[];
low_bou=[-2;-2;-2;-2];
up_bou=[2;2;2;2];
nonlcon_function=[];
cheapcon_function=[];

% variable_number=2;
% object_function_H=@(x) func_H2_2D_H(x);
% object_function_L=@(x) func_H2_2D_L(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[-5;-5];
% up_bou=[5;5];
% nonlcon_function=[];
% cheapcon_function=[];

x_intial=rand(variable_number,1).*(up_bou-low_bou)+low_bou;
[x_best_FM,fval_best_FM]=fmincon(object_function_H,x_intial,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function)

delete([data_library_name,'.txt']);
delete([data_library_name,'_source.txt']);
delete('result_total.txt');

model_function=@(x) modelFunction(x,object_function_L,nonlcon_function);
[~,x_list,~]=getLatinHypercube...
    (min(100,10*variable_number),variable_number,[],...
    low_bou,up_bou,cheapcon_function);
dataLibraryUpdata...
    ([data_library_name,'_source'],model_function,x_list);

[x_best,fval_best,NFE,output]=optimalSurrogateSADEKTS...
    (object_function_H,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,[])
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
    (object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,NFE_max)
% surrogate base optimal use polynomial responed surface method version 1
% use MSP guideline to decide which point was select to updata, ga optimal
% iteration the model an optimal
% x_best is row vector by ga, so need to transpose into colume vector
% all function exchange x should be colume vector,x_list except
% both nonlcon_function and cheapcon_function format is [con,coneq]
% model_function should output fval, format is [fval,con,coneq]
% con or coneq can be colume vector if there was more than one constrain
%
% Copyright 2022 Adel
%
if nargin < 12
    NFE_max=50*variable_number;
    if nargin <10
        cheapcon_function=[];
        if nargin < 9
            nonlcon_function=[];
        end
    end
end
INFORMATION_FLAG=0; % whether draw data
CONVERGENCE_JUDGMENT_FLAG=0; % whether judgment convergence

population_number=min(100,10*variable_number);
correction_factor=0.1; % ahpla
elite_rate=0.5;

scaling_factor=0.8; % F
cross_rate=0.8;

iteration_max=100;
torlance=1e-4;

data_library_name='optimalSurrogateSADEKTS_result';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;
result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);

% if do not input model_function, generate model_function
if ~isempty(nonlcon_function)
    expensive_nonlcon_flag=1;
else
    expensive_nonlcon_flag=0;
end
if nargin < 11 || isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end
if isempty(A) && isempty(B) && isempty(Aeq) && isempty(Beq) &&...
        isempty(cheapcon_function)
    cheapcon_function=[];
else
    cheapcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq,cheapcon_function);
end

iteration=iteration+1;

% step 2
% generate initial sample x_list
X_updata=getInitialSample...
    (population_number,[data_library_name,'_source'],variable_number,low_bou,up_bou,elite_rate,correction_factor);

while ~done
    % step 3
    % updata data library by x_list
    dataLibraryUpdata(data_library_name,model_function,X_updata);NFE=NFE+size(X_updata,1);
    
    % import data from data library to rank data
    [x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
        (data_library_name,low_bou,up_bou);
    
    % nomalization con
    fval_max=max(abs(fval_list));
    if ~isempty(con_list)
        con_max_list=max(abs(con_list),[],1);
        con_list=con_list./con_max_list.*fval_max;
    end
    if ~isempty(coneq_list)
        coneq_max_list=max(abs(coneq_list),[],1);
        coneq_list=coneq_list./coneq_max_list.*fval_max;
    end
    
    % step 4
    % rank x_list data
    [x_list,fval_list,con_list,coneq_list]=rankData...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function);
    
    % step 5
    % only the first population_number will be use
    x_list_best_po=x_list(1:population_number,:);
    
    % differ evolution mutations
    X_new_R1=differEvolutionRand...
        (low_bou,up_bou,x_list_best_po,scaling_factor,population_number,1);
    X_new_R2=differEvolutionRand...
        (low_bou,up_bou,x_list_best_po,scaling_factor,population_number,2);
    X_new_CR=differEvolutionCurrentRand...
        (low_bou,up_bou,x_list_best_po,scaling_factor);
    X_new_CB=differEvolutionCurrentBest...
        (low_bou,up_bou,x_list_best_po,scaling_factor,1);
    
    % differ evolution crossover
    X_new_R1=differEvolutionCrossover...
        (low_bou,up_bou,x_list_best_po,X_new_R1,cross_rate);
    X_new_R2=differEvolutionCrossover...
        (low_bou,up_bou,x_list_best_po,X_new_R2,cross_rate);
    X_new_CR=differEvolutionCrossover...
        (low_bou,up_bou,x_list_best_po,X_new_CR,cross_rate);
    X_new_CB=differEvolutionCrossover...
        (low_bou,up_bou,x_list_best_po,X_new_CB,cross_rate);
    
    % step 6
    % get kriging model and function
    [kriging_model_fval,~,~,output_Kriging]=getKrigingModel...
        (x_list,fval_list,con_list,coneq_list);
    object_function=output_Kriging.object_function;
    nonlcon_function=output_Kriging.nonlcon_function;
    
    % find global infill point base kriging model from offspring X
    x_list_offspring=[X_new_R1;X_new_R2;X_new_CR;X_new_CB];
    fval_list_offspring=zeros(4*population_number,1);
    con_list_offspring=zeros(4*population_number,size(con_list,2));
    coneq_list_offspring=zeros(4*population_number,size(coneq_list,2));
    for x_index=1:4*population_number
        x=x_list_offspring(x_index,:)';
        fval_list_offspring(x_index,1)=object_function(x);
        if ~isempty(nonlcon_function)
            [con,coneq]=nonlcon_function(x);
            if ~isempty(con)
                con_list_offspring(x_index,:)=con';
            end
            if ~isempty(coneq)
                coneq_list_offspring(x_index,:)=coneq';
            end
        end
    end
    [x_list_offspring,~,~,~]=rankData...
        (x_list_offspring,fval_list_offspring,con_list_offspring,coneq_list_offspring,...
        cheapcon_function);
    x_global_infill=x_list_offspring(1,:)';
    
    % distance to exist point can not be to low
    distance=min(sum((x_list-x_global_infill').^2./(low_bou'-up_bou').^2,2));
    if  distance > torlance
        % updata global infill point into data library
        [~,~,~]=dataLibraryUpdata(data_library_name,model_function,x_global_infill');
    end
    
    % step 7
    % rand select local point from X_g
    x_index=randi(4*population_number);
    x_initial=x_list_offspring(x_index,:)';
    %     [x_best,~]=findMinRaw...
    %         (x_list,fval_list,con_list,coneq_list,...
    %         A,B,Aeq,Beq,cheapcon_function);
    %     x_initial=x_best;
    
    % select nearest point to construct RBF
    distance=sum(x_initial'-x_list,2);
    [~,index]=sort(distance);
    if (variable_number+1)*(variable_number+2)/2 >= 100
        RBF_number=min((variable_number+1)*(variable_number+2)/2,size(x_list,1));
    else
        RBF_number=min(100,size(x_list,1));
    end
    index=index(1:RBF_number);
    x_list_RBF=x_list(index,:);
    fval_list_RBF=fval_list(index,:);
    if ~isempty(con_list)
        con_list_RBF=con_list(index,:);
    else
        con_list_RBF=[];
    end
    if ~isempty(coneq_list)
        coneq_list_RBF=coneq_list(index,:);
    else
        coneq_list_RBF=[];
    end
    
    % get RBF model and function
    [radialbasis_model_fval,~,~,output_radialbasis]=getRadialBasisModel...
        (x_list_RBF,fval_list_RBF,con_list_RBF,coneq_list_RBF);
    object_function=output_radialbasis.object_function;
    nonlcon_function=output_radialbasis.nonlcon_function;
    
    if INFORMATION_FLAG && variable_number < 3
        figure(5);
        scatter(x_list(:,1),x_list(:,2));
        line(x_global_infill(1),x_global_infill(2),'Marker','o','color','r')
        interpolationKrigingVisualize(kriging_model_fval,low_bou,up_bou);
        interpolationRadialBasisVisualize(radialbasis_model_fval,low_bou,up_bou);
    end
    
    % get local infill point
    fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp');
    x_local_infill=fmincon(object_function,x_initial,A,B,Aeq,Beq,...
        low_bou,up_bou,nonlcon_function,fmincon_options);
    
    % recond min fval each iteration
    [x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
        (data_library_name,low_bou,up_bou);
    [x_best,fval_best,con_best,coneq_best]=findMinRaw...
        (x_list,fval_list,con_list,coneq_list,...
        cheapcon_function);
    result_x_best(iteration,:)=x_best';
    result_fval_best(iteration,:)=fval_best;
    iteration=iteration+1;
    
    % forced interrupt
    if iteration>=iteration_max || NFE >=NFE_max
        done=1;
    end
    
    % convergence judgment
    if CONVERGENCE_JUDGMENT_FLAG
        if (iteration > iteration_max*0.2 && ...
                abs((fval_potential-fval_potential_old)/fval_potential_old) < torlance)
            if ~sum(con_best > torlance)
                done=1;
            end
        end
    end
    
    % distance to exist point can not be to low
    distance=min(sum((x_list-x_local_infill').^2./(low_bou'-up_bou').^2,2));
    if  distance > torlance
        % updata global infill point into data library
        X_updata=x_local_infill';
    else
        X_updata=[];
    end
    
end

result_x_best=result_x_best(1:iteration-1,:);
result_fval_best=result_fval_best(1:iteration-1);

output.result_x_best=result_x_best;
output.result_fval_best=result_fval_best;

    function [con,coneq]=cheapconFunction(x,A,B,Aeq,Beq,cheapcon_function)
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
                index_bou__=find(X_new(x_index__,:) < low_bou');
                X_new(x_index__,index_bou__)=low_bou(index_bou__)';
                index_bou__=find(X_new(x_index__,:) > up_bou');
                X_new(x_index__,index_bou__)=up_bou(index_bou__)';
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
            index_bou__=find(X_new(x_index__,:) < low_bou');
            X_new(x_index__,index_bou__)=low_bou(index_bou__)';
            index_bou__=find(X_new(x_index__,:) > up_bou');
            X_new(x_index__,index_bou__)=up_bou(index_bou__)';
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
            index_bou__=find(X_new(x_index__,:) < low_bou');
            X_new(x_index__,index_bou__)=low_bou(index_bou__)';
            index_bou__=find(X_new(x_index__,:) > up_bou');
            X_new(x_index__,index_bou__)=up_bou(index_bou__)';
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
        index_bou__=find(X_new < low_bou');
        X_new(index_bou__)=low_bou(index_bou__)';
        index_bou__=find(X_new > up_bou');
        X_new(index_bou__)=up_bou(index_bou__)';
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

function [x_best,fval_best,con_best,coneq_best]=findMinRaw...
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
con_best=[];
coneq_best=[];
if ~isempty(con_list)
    con_sum_list=con_sum_list+sum(max(con_list,0),2);
end
if ~isempty(coneq_list)
    con_sum_list=con_sum_list+sum(coneq_list.^2,2);
end

% add cheap con
for x_index=1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x_list(x_index,:)');
        con_sum_list(x_index)=con_sum_list(x_index)+...
            sum(max(con,0))+sum(coneq.*coneq);
    end
end

index=find(con_sum_list<=0);
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
    [~,index_best]=min(con_sum_list);
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

function [x_list]=getInitialSample...
    (population_number,data_library_name,variable_number,...
    low_bou,up_bou,...
    elite_rate,correction_factor)
% use SVM to correct latin hypercubic by exist data
% Knowledge-Transfer-Based Sampling Method
%

% generate initial latin hypercubic
[X,X_new,distance_min_normalize]=getLatinHypercube...
    (population_number,variable_number,[],...
    low_bou,up_bou);

% import data from data library to rank data
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    ([data_library_name,'.txt'],low_bou,up_bou);

% KTS
if ~isempty(x_list)
    % rank x_list data
    [X_input,~,~,~]=rankData...
        (x_list,fval_list,con_list,coneq_list);
    
    [x_number,variable_number]=size(x_list);
    
    N_elite=round(x_number*elite_rate);
    
    % generate SVM model
    Y=[ones(N_elite,1);-ones(x_number-N_elite,1)];
    kernel_function=@(x1,x2) exp(-((x1-x2)'*(x1-x2))/2*1000);
    [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
        (X_input,Y,10,kernel_function);
    
    % get predict value by SVM
    Y=zeros(population_number,1);
    index_list=[];
    for x_index=1:population_number
        Y(x_index)=SVM_predict_function(X(x_index,:)');
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
    X_superior=X(index_list,:);
    X(index_list,:)=[];
    x_list=[];
    for x_index=1:size(X,1)
        x=X(x_index,:);
        distance=sqrt(sum((x-X_superior).^2,2));
        [~,index]=min(distance);
        x_superior=X_superior(index,:); % nearest x_superior
        x=x+correction_factor*(x_superior-x);
        x_list=[x_list;x];
    end
    
    x_list=[x_list;X_superior];
else
    x_list=X;
end
end
function [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
    (X,Y,C,kernel_function,low_bou,up_bou)
% generate support vector machine model
% X, Y, colume is each iterm
% low_bou, up_bou is colume vector
% C is penalty factor, default is empty
% kernel_function default is X*XT
% astdD_X, stdD_X is row vector
% kernel_function input colume vector
%
if nargin < 4
    kernel_function=[];
    if nargin < 3
        C=[];
    end
end

[x_number,variable_number]=size(X);

astdD_X=mean(X);
stdD_X=std(X);
% normalization data
if nargin < 6
    up_bou=max(X)';
    if nargin < 5
        low_bou=min(X)';
    end
end
X_norml=(X-low_bou')./(up_bou'-low_bou');

if isempty(kernel_function)
    kernel_function=@(x1,x2) x1'*x2;
end

X_inner_product=zeros(x_number);
for index_i=1:x_number
    for index_j=1:x_number
        X_inner_product(index_i,index_j)=...
            kernel_function(X_norml(index_i,:)',X_norml(index_j,:)');
    end
end

% min SVM object function to get lamada
object_function_SVM=@(lamada) -objectFunctionSVM(lamada,X_inner_product,Y);
lamada_initial=ones(x_number,1)*0.5;
low_bou_fmincon=0*ones(x_number,1);
if isempty(C) || C==0
    up_bou_fmincon=[];
else
    up_bou_fmincon=C*ones(x_number,1);
end
Aeq=Y';
fmincon_options = optimoptions('fmincon','Display','none','Algorithm','sqp');

lamada=fmincon(object_function_SVM,lamada_initial,...
    [],[],Aeq,0,low_bou_fmincon,[],[],fmincon_options);

% obtain other paramter
w=sum(lamada.*Y.*X_norml);
index_list=find(lamada > 1e-6);
b=sum(Y(index_list)-...
    (sum(lamada.*Y.*X_inner_product(:,index_list),1))')...
    /size(index_list,1);

% disp('error');

% generate predict function
SVM_predict_function=@(x) classifySupportVectorMachinePredictor...
    (x,X_norml,Y,lamada,b,low_bou,up_bou,kernel_function);

% output model
SVM_model.X=X;
SVM_model.Y=Y;
SVM_model.X_norml=X_norml;
SVM_model.low_bou=low_bou;
SVM_model.up_bou=up_bou;
SVM_model.astdD_X=astdD_X;
SVM_model.stdD_X=stdD_X;
SVM_model.lamada=lamada;
SVM_model.w=w;
SVM_model.b=b;
SVM_model.kernel_function=kernel_function;
SVM_model.SVM_predict_function=SVM_predict_function;

    function fval=objectFunctionSVM(lamada,X_inner_product,Y)
        % support vector machine maximum object function
        %
        lamada=lamada(:);
        lamada_Y=lamada.*Y;
        fval=sum(lamada)-lamada_Y'*X_inner_product*lamada_Y/2;
    end

    function [y,y_sum]=classifySupportVectorMachinePredictor...
            (x,X_norml,Y,lamada,b,low_bou,up_bou,kernel_function)
        % predict value of x is 1 or -1
        % x input is colume vector
        %
        x_number__=size(X_norml,1);
        x_norml=(x-low_bou)./(up_bou-low_bou);
        X_inner_product__=zeros(x_number__,1);
        for x_index__=1:x_number__
            X_inner_product__(x_index__)=...
                kernel_function(X_norml(x_index__,:)',x_norml);
        end
        y_sum=sum(lamada.*Y.*X_inner_product__)+b;
        if y_sum > 0
            y=1;
        elseif y_sum == 0
            y=0;
        else
            y=-1;
        end
    end
end
function classifySupportVectorMachineVisualization...
    (SVM_model,low_bou,up_bou)
% Visualization SVM_model
%
if nargin < 1
    error('classifySupportVectorMachineVisualization: not enough input');
end
X=SVM_model.X;
Y=SVM_model.Y;
% norml_X=SVM_model.norml_X;
% lamada=SVM_model.lamada;
% w=SVM_model.w;
% b=SVM_model.b;
% kernel_function=SVM_model.kernel_function;
SVM_predict_function=SVM_model.SVM_predict_function;

if nargin < 3
    up_bou=max(X)';
    if nargin < 2
        low_bou=min(X)';
    end
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
[X,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
fval=zeros(grid_number+1);
fval_sum=zeros(grid_number+1);
for x_index=1:grid_number+1
    for y_index=1:grid_number+1
        predict_x=([x_index;y_index]-1).*d_bou+low_bou;
        [fval(y_index,x_index),fval_sum(y_index,x_index)]=...
            SVM_predict_function(predict_x);
    end
end
contour(X,Y,fval,[0,0]);

% draw point
line(X_positive(:,1),X_positive(:,2),'LineStyle','none','Marker','o','Color','b');
line(X_negative(:,1),X_negative(:,2),'LineStyle','none','Marker','o','Color','r');

% figure(2)
% surf(X,Y,fval_sum);
end

function [x_list,fval_list,con_list,coneq_list]=rankData...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function)
% rank data base on feasibility rule
% infeasible is rank by sum of constraint
%
if nargin < 9
    cheapcon_function=[];
    if nargin < 8
        Beq=[];
        if nargin < 7
            Aeq=[];
            if nargin < 6
                B=[];
                if nargin < 5
                    A=[];
                end
            end
        end
    end
end

[x_number,~]=size(x_list);
con_sum_list=zeros(x_number,1);
for x_index=1:x_number
    con_sum_list(x_index)=sum(max(con_list(x_index,:),0))+...
        sum(abs(coneq_list(x_index,:)));
end

% add cheap con
for x_index=1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq]=cheapcon_function(x_list(x_index,:)');
        con_sum_list(x_index)=con_sum_list(x_index)+...
            sum(max(con,0))+sum(max(coneq,0));
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
    (x_list,fval_list,con_list,coneq_list,...
    kriging_model_fval,kriging_model_con,kriging_model_coneq)
% base on library_data to create kriging model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
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
% base function£¬which can be a p x one function matrix
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
    figure_handle=figure(103);
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

function [radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq,output]=getRadialBasisModel...
    (x_list,fval_list,con_list,coneq_list,...
    radialbasis_model_fval,radialbasis_model_con,radialbasis_model_coneq)
% base on library_data to create radialbasis model and function
% if input model, function will updata model
% object_function is single fval output
% nonlcon_function is normal nonlcon_function which include con, coneq
% con is colume vector, coneq is colume vector
% var_function is same
%
if nargin < 7
    radialbasis_model_coneq=[];
    if nargin < 6
        radialbasis_model_con=[];
        if nargin < 5
            radialbasis_model_fval=[];
        end
    end
end

% generate surrogate model
if isempty(radialbasis_model_fval)
    radialbasis_model_fval=interpolationRadialBasisPreModel...
        (x_list,fval_list);
else
    radialbasis_model_fval=interpolationRadialBasisUpdata...
        (radialbasis_model_fval,x_list,fval_list);
end

if ~isempty(con_list)
    if isempty(radialbasis_model_fval)
        radialbasis_model_con=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
            'basis_function',[],'beta',[],...
            'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[]);
        radialbasis_model_con=repmat(radialbasis_model_con,[size(con_list,2),1]);
        for con_index=1:size(con_list,2)
            radialbasis_model_con(con_index)=interpolationRadialBasisPreModel...
                (x_list,con_list(:,con_index));
        end
    else
        for con_index=1:size(con_list,2)
            radialbasis_model_con(con_index)=interpolationRadialBasisUpdata...
                (radialbasis_model_con(con_index),x_list,con_list(:,con_index));
        end
    end
end

if ~isempty(coneq_list)
    if isempty(radialbasis_model_coneq)
        radialbasis_model_coneq=struct('X',[],'Y',[],'X_normalize',[],'Y_normalize',[],...
            'basis_function',[],'beta',[],...
            'aver_X',[],'stdD_X',[],'aver_Y',[],'stdD_Y',[]);
        radialbasis_model_coneq=repmat(radialbasis_model_coneq,[size(coneq_list,2),1]);
        for coneq_index=1:size(coneq_list,2)
            radialbasis_model_coneq(coneq_index)=interpolationRadialBasisPreModel...
                (x_list,coneq_list(:,coneq_index));
        end
    else
        for coneq_index=1:size(coneq_list,2)
            radialbasis_model_coneq(coneq_index)=interpolationRadialBasisUpdata...
                (radialbasis_model_coneq(coneq_index),x_new_list,coneq_list(:,coneq_index));
        end
    end
end

object_function=@(predict_x) objectFunctionSurrogate(predict_x,radialbasis_model_fval);
if isempty(radialbasis_model_con) && isempty(radialbasis_model_coneq)
    nonlcon_function=[];
else
    nonlcon_function=@(predict_x) constraintFunctionSurrogate(predict_x,radialbasis_model_con,radialbasis_model_coneq);
end

output.object_function=object_function;
output.nonlcon_function=nonlcon_function;
output.x_list=x_list;
output.fval_list=fval_list;
output.con_list=con_list;
output.coneq_list=coneq_list;


    function fval=objectFunctionSurrogate...
            (predict_x,radialbasis_model_fval)
        fval=interpolationRadialBasisPredictor(radialbasis_model_fval,predict_x);
    end
    function [con,coneq]=constraintFunctionSurrogate...
            (predict_x,radialbasis_model_con,radialbasis_model_coneq)
        if isempty(radialbasis_model_con)
            con=[];
        else
            con=zeros(length(radialbasis_model_con),1);
            for con_ind=1:length(radialbasis_model_con)
                con(con_ind)=interpolationRadialBasisPredictor...
                    (radialbasis_model_con(con_ind),predict_x);
            end
        end
        if isempty(radialbasis_model_coneq)
            coneq=[];
        else
            coneq=zeros(length(radialbasis_model_coneq),1);
            for coneq_ind=1:length(radialbasis_model_con)
                coneq(coneq_ind)=interpolationRadialBasisPredictor...
                    (radialbasis_model_coneq(coneq_ind),predict_x);
            end
        end
    end
end
function radialbasis_model=interpolationRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interpolation pre model function
% input initial data X and Y, x is row vector
% return model
% beta in radialbasis_model is colume vector
%
% Copyright 2022 Adel
%
if nargin < 3
    basis_function=[];
end
if isempty(basis_function)
    basis_function=@(r) exp(-(r'*r)/2);
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
X_norml = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_norml = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

beta=interpolationRadialBasis(X_norml,Y_norml,basis_function);

radialbasis_model.X=X;
radialbasis_model.Y=Y;
radialbasis_model.X_normalize=X_norml;
radialbasis_model.Y_normalize=Y_norml;
radialbasis_model.beta=beta;
radialbasis_model.aver_X=aver_X';
radialbasis_model.stdD_X=stdD_X';
radialbasis_model.aver_Y=aver_Y';
radialbasis_model.stdD_Y=stdD_Y';
radialbasis_model.basis_function=basis_function;
end
function [predict_y]=interpolationRadialBasisPredictor...
    (radialbasis_model,predict_x)
% radial basis function interpolation predict function
% input predict_x and radialbasis_model model
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
X_norml=radialbasis_model.X_normalize;
beta=radialbasis_model.beta;
aver_X=radialbasis_model.aver_X;
stdD_X=radialbasis_model.stdD_X;
aver_Y=radialbasis_model.aver_Y;
stdD_Y=radialbasis_model.stdD_Y;
basis_function=radialbasis_model.basis_function;

[x_number,variable_number]=size(X_norml);

% normalize data
predict_x=(predict_x-aver_X')./stdD_X';

% predict value
X_inner_product=zeros(x_number,1);
for index_i=1:x_number
    X_inner_product(index_i,1)=...
        basis_function(X_norml(index_i,:)'-predict_x');
end

% predict variance
predict_y=beta'*X_inner_product;

% normalize data
predict_y=predict_y*stdD_Y+aver_Y;
end
function radialbasis_model=interpolationRadialBasisUpdata...
    (radialbasis_model,X_new,Y_new)
% add new point in existing radial basis function model
% X_new is row vector, can be more than one
% Y_new is row vector, can be more than one
%

% giving value
X=radialbasis_model.X;
Y=radialbasis_model.Y;
basis_function=radialbasis_model.basis_function;

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
beta=interpolationRadialBasis(X_normalize,Y_normalize,basis_function);

radialbasis_model.X=X;
radialbasis_model.Y=Y;
radialbasis_model.X_normalize=X_normalize;
radialbasis_model.Y_normalize=Y_normalize;
radialbasis_model.beta=beta;
radialbasis_model.aver_X=aver_X';
radialbasis_model.stdD_X=stdD_X';
radialbasis_model.aver_Y=aver_Y';
radialbasis_model.stdD_Y=stdD_Y';
end
function interpolationRadialBasisVisualize...
    (radialbasis_model,low_bou,up_bou,figure_handle)
% visualization polynamial respond surface model
% figrue is 100
%
if nargin < 4
    figure_handle=figure(102);
end
if size(low_bou,1) ~= size(low_bou,1)
    error('interpolationRadialBasisVisualize: boundary incorrect');
end
if size(low_bou,1) > 2
    error('interpolationRadialBasisVisualize: dimension large than two');
end
axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end

x_list=radialbasis_model.X;
y_list=radialbasis_model.Y;

grid_number=100;
d_bou=(up_bou-low_bou)/grid_number;
predict_result=zeros(grid_number+1);
[X,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
for x_index=1:grid_number+1
    for y_index=1:grid_number+1
        predict_x=([x_index,y_index]-1).*d_bou'+low_bou';
        [predict_result(y_index,x_index)]=interpolationRadialBasisPredictor...
            (radialbasis_model,predict_x);
    end
end
surf(axes_handle,X,Y,predict_result,'FaceAlpha',0.5,'EdgeColor','none');
line(axes_handle,x_list(:,1),x_list(:,2),y_list,'Marker','o','LineStyle','none');
xlabel('X');
ylabel('Y');
zlabel('Z');
end
function beta=interpolationRadialBasis...
    (X,Y,basis_function)
% interpolation polynomial responed surface core function
% calculation beta
%

[x_number,variable_number]=size(X);

X_inner_product=zeros(x_number);
for index_i=1:x_number
    for index_j=1:x_number
        X_inner_product(index_i,index_j)=...
            basis_function(X(index_i,:)'-X(index_j,:)');
    end
end

beta=X_inner_product\Y;
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
file_optimalSurrogate_output = fopen([data_library_name,'.txt'],'a');
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

if exist([data_library_name,'.txt'],'file')==2
    data_list=importdata([data_library_name,'.txt']);
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
    cheapcon_function=@(x) ...
        max(max(sample_number*cheapcon_function(x.*(up_bou-low_bou)+low_bou)+1,0),[],1);
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
                    fval__(x_index)=object_function(X__(x_index,:));
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
                if max(cheapcon_function(x_supply),[],1) > 0
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