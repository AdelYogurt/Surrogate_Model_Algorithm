clc;
clear;
close all hidden;

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

% variable_number=2;
% object_function=@functionG06Object;
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[13;0];
% up_bou=[100;100];
% cheapcon_function=[];
% nonlcon_function=@functionG06Nonlcon;
% model_function=[];

% variable_number=13;
% object_function=@func_G3_G01;
% A=[
%     2   2   0   0   0   0   0   0   0   1   1   0   0;
%     2   0   2   0   0   0   0   0   0   1   0   1  0;
%     0   2   2   0   0   0   0   0   0   0   1   1  0;
%     -8  0   0   0   0   0   0   0   0   1   0   0   0;
%     0   -8  0   0   0   0   0   0   0   0   1   0   0;
%     0   0   0   -2  -1  0   0   0   0   1   0   0   0;
%     0   0   0   0   0   -2  -1  0   0   0   1   0   0;
%     0   0   0   0   0   0   0   -2  -1  0   0   1   0;
% ];
% B=[10 10 10 0 0 0 0 0]';
% Aeq=[];
% Beq=[];
% low_bou=zeros(13,1);
% up_bou=ones(13,1);
% up_bou(10:12)=100;
% nonlcon_function=[];
% cheapcon_function=[];

% variable_number=4;
% object_function=@(x) functionPVDObject(x);
% object_function_low=@(x) functionPVDObjectLow(x);
% A=[-1,0,0.0193,0;
%     0,-1,0.00954,0;];
% B=[0;0];
% Aeq=[];
% Beq=[];
% low_bou=[1.0;0.625;25;25];
% up_bou=[1.375;1;150;240];
% nonlcon_function=@(x) functionPVDNonlcon(x);
% cheapcon_function=@(x) cheapconFunction(x,A,B,Aeq,Beq);
% model_function=[];

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

% variable_number=10;
% object_function=@(x) functionR10Object(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=zeros(variable_number,1);
% up_bou=ones(variable_number,1);
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];

% variable_number=5;
% object_function=@(x) functionICEObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=[50;2;20;20;4000];
% up_bou=[1e2;10;50;50;10000];
% nonlcon_function=@(x) functionICENonlcon(x);
% cheapcon_function=[];

% nonlcon_function=[];
% cheapcon_function=@nonlcon_ICE;

% x_initial=rand(variable_number,1).*(up_bou-low_bou)+low_bou;
% [x_best,fval_best]=fmincon(object_function,x_initial,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function)


data_library_name='optimalSurrogate_PRS_TR_result.txt';
delete(data_library_name);
delete('result_total.txt');
[x_best,fval_best,NFE,output]=optimalSurrogatePRSTR...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,[],100)
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


% repeat_time=50;
% fval_best_list=zeros(repeat_time,1);
% for repeat_index=1:repeat_time
%     delete(data_library_name);
%     delete('result_total.txt');
%     [x_best,fval_best,NFE,output]=optimalSurrogateRespSurf...
%         (object_function,variable_number,A,B,Aeq,Beq,low_bou,up_bou,nonlcon_function,...
%         cheapcon_function,[]);
%     fval_best_list(repeat_index,1)=fval_best;
% end
% boxplot(fval_best_list);

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

function [x_best,fval_best,NFE,output]=optimalSurrogatePRSTR...
    (object_function,variable_number,low_bou,up_bou,nonlcon_function,...
    cheapcon_function,model_function,....
    NFE_max,iteration_max,torlance,nonlcon_torlance)
% surrogate base optimal  Trust-Region-Based Adaptive Response Surface version 3
% use MSP guideline to decide which point was select to updata, spq optimal
% x_best is row vector by ga, so need to transpose into colume vector
% all function exchange x should be colume vector,x_list except
% both nonlcon_function and cheapcon_function format is [con,coneq]
% model_function should output fval, format is [fval,con,coneq]
% con or coneq can be colume vector if there was more than one constrain
%
% Copyright 2022 10 Adel
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

% if noob user input cheapcon function handle while no constraint
if ~isempty(cheapcon_function)
    [con,coneq]=cheapcon_function(rand(variable_number,1).*(up_bou-low_bou)+low_bou);
    if isempty(con) && isempty(coneq)
        cheapcon_function=[];
    end
end
DRAW_FIGURE_FLAG=0; % whether draw data
INFORMATION_FLAG=0; % whether print data
CONVERGENCE_JUDGMENT_FLAG=0; % whether judgment convergence

% parameter
enlarge_range=2; % adapt region enlarge parameter
range_max=0.5;

protect_range=1e-4; % point added into data library min range

% augmented lagrange parameter
lambda_initial=10;
miu=1;
miu_max=1000;
gama=2;

sample_number=(variable_number+1)*(variable_number+2)/2; % Latin hypercube sample count

data_library_name='optimalSurrogate_PRS_TR_result.txt';
file_result = fopen('result_total.txt','a');
fprintf(file_result,'%s\n',datetime);
fclose(file_result);
clear('file_result');

done=0;NFE=0;iteration=0;

% if do not input model_function, generate model_function
if nargin < 7 || isempty(model_function)
    model_function=@(x) modelFunction(x,object_function,nonlcon_function);
end

iteration=iteration+1;

% bourdary updata
low_bou_iter=low_bou;
up_bou_iter=up_bou;

% Step 1
% generate latin hypercube sequence
[x_list_exist,~,con,coneq]=dataLibraryLoad...
    (data_library_name,low_bou_iter,up_bou_iter);
[~,x_updata_list,~]=getLatinHypercube...
    (sample_number,variable_number,x_list_exist,...
    low_bou_iter,up_bou_iter,cheapcon_function);

% if only input model function, detect constraint
if ~isempty(x_updata_list)
    [~,con,coneq]=dataLibraryUpdata...
        (data_library_name,model_function,x_updata_list(1,:));NFE=NFE+1;
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
x_updata_list=x_updata_list(2:end,:);
result_x_best=zeros(iteration_max,variable_number);
result_fval_best=zeros(iteration_max,1);

% total data storage in x_list, fval_list, con_list, coneq_list
[x_list,fval_list,con_list,coneq_list]=dataLibraryLoad...
    (data_library_name,low_bou_iter,up_bou_iter);

% loop
while ~done
    % Step 2
    % update x_updata_list into data library
    [fval_updata_list,con_updata_list,coneq_updata_list]=dataLibraryUpdata...
        (data_library_name,model_function,x_updata_list);NFE=NFE+size(x_updata_list,1);
    x_list=[x_list;x_updata_list];
    fval_list=[fval_list;fval_updata_list];
    if ~isempty(con_list)
        con_list=[con_list;con_updata_list];
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_updata_list];
    end
    
    % Step 3
    % load data
    % RPS_list is used to interpolation
    [x_RPS_list,fval_RPS_list,con_RPS_list,coneq_RPS_list]=dataLibraryLoad...
        (data_library_name,low_bou_iter,up_bou_iter);
    
    %     % nomalization con
    %     fval_max=max(abs(fval_list));
    %     if ~isempty(con_list)
    %         con_max_list=max(abs(con_list),[],1);
    %         con_list=con_list./con_max_list.*fval_max;
    %     end
    %     if ~isempty(coneq_list)
    %         coneq_max_list=max(abs(coneq_list),[],1);
    %         coneq_list=coneq_list./coneq_max_list.*fval_max;
    %     end
    
    % nomalization con
    fval_RPS_max=mean(abs(fval_RPS_list),1);
    fval_RPS_list=fval_RPS_list./fval_RPS_max*1e3;
    if ~isempty(con_RPS_list)
        con_RPS_max_list=mean(abs(con_RPS_list),1);
        con_RPS_list=con_RPS_list./con_RPS_max_list*1e3;
    end
    if ~isempty(coneq_RPS_list)
        coneq_RPS_max_list=mean(abs(coneq_RPS_list),1);
        coneq_RPS_list=coneq_RPS_list./coneq_RPS_max_list*1e3;
    end
    
    % get surrogate model and function
    % avoid too close point
    select_index_list=1:size(x_RPS_list,1);
    index=1;
    while index < length(select_index_list)
        select_index=select_index_list(index);
        distance=sum(((x_RPS_list(select_index,:)-...
            [x_RPS_list(1:select_index-1,:);x_RPS_list(select_index+1:end,:)])./...
            (up_bou_iter'-low_bou_iter')).^2,2);
        if min(distance) < protect_range^2
            select_index_list(index)=[];
        else
            index=index+1;
        end
    end
    x_RPS_list=x_RPS_list(select_index_list,:);
    fval_RPS_list=fval_RPS_list(select_index_list,:);
    if ~isempty(con_RPS_list)
        con_RPS_list=con_RPS_list(select_index_list,:);
    end
    if ~isempty(coneq_RPS_list)
        coneq_RPS_list=coneq_RPS_list(select_index_list,:);
    end
    
    % if point too less, add more point
    if size(x_RPS_list,1) < (variable_number+1)*(variable_number+2)/2
        % generate latin hypercube sequence
        [~,x_updata_list,~]=getLatinHypercube...
            (sample_number,variable_number,x_RPS_list,...
            low_bou_iter,up_bou_iter,cheapcon_function);
        % update x_updata_list into data library
        [fval_updata_list,con_updata_list,coneq_updata_list]=dataLibraryUpdata...
            (data_library_name,model_function,x_updata_list);NFE=NFE+size(x_updata_list,1);
        x_list=[x_list;x_updata_list];
        fval_list=[fval_list;fval_updata_list];
        con_list=[con_list;con_updata_list];
        coneq_list=[coneq_list;coneq_updata_list];
        
        % normalization data and updata into list
        x_RPS_list=[x_RPS_list;x_updata_list];
        fval_updata_nomlz_list=fval_updata_list./fval_RPS_max*1e3;
        fval_RPS_list=[fval_RPS_list;fval_updata_nomlz_list];
        if ~isempty(con_RPS_list)
            con_updata_list=con_updata_list./con_RPS_max_list*1e3;
            con_RPS_list=[con_RPS_list;con_updata_list];
        end
        if ~isempty(coneq_RPS_list)
            coneq_updata_list=coneq_updata_list./coneq_RPS_max_list*1e3;
            coneq_RPS_list=[coneq_RPS_list;coneq_updata_list];
        end
        
    end
    [respsurf_model_fval,~,~,output_respsurf]=getRespSurfModel...
        (x_RPS_list,fval_RPS_list,con_RPS_list,coneq_RPS_list);
    object_function_surrogate=output_respsurf.object_function_surrogate;
    nonlcon_function_surrogate=output_respsurf.nonlcon_function_surrogate;
    
    % generate merit function
    if expensive_nonlcon_flag
        if iteration == 1
            con_number=size(con_list,2);
            coneq_number=size(coneq_list,2);
            lambda_con=lambda_initial*ones(con_number,1);
            lambda_coneq=lambda_initial*ones(coneq_number,1);
        end
        % generate merit function fval and add fval
        merit_function=@(x) meritFunction...
            (x,object_function_surrogate,nonlcon_function_surrogate,...
            lambda_con,lambda_coneq,miu);
    end
    
    if DRAW_FIGURE_FLAG
        interpolationVisualize(respsurf_model_fval,low_bou_iter,up_bou_iter);
%         drawFunction(merit_function,low_bou_iter,up_bou_iter)
    end
    
    % Step 4
    % updata data, x_potential
    [x_best,~,~,~]=findMinRaw...
        (x_RPS_list,fval_RPS_list,con_RPS_list,coneq_RPS_list,...
        cheapcon_function,nonlcon_torlance);
    if expensive_nonlcon_flag
        [x_potential,fval_potential_predict]=findMinMSPFmincon...
            (merit_function,x_best,low_bou_iter,up_bou_iter,...
            cheapcon_function);
    else
        [x_potential,fval_potential_predict]=findMinMSPFmincon...
            (object_function_surrogate,x_best,low_bou_iter,up_bou_iter,...
            cheapcon_function);
    end
    
    % check x_potential if exist in data library
    % if not, updata data libraray
    [x_potential,fval_potential,con_potential,coneq_potential,NFE_prot]=dataLibraryUpdataProtect...
        (data_library_name,model_function,x_potential',...
        x_list,fval_list,con_list,coneq_list,...
        low_bou,up_bou,protect_range);NFE=NFE+NFE_prot;
    x_list=[x_list;x_potential];x_potential=x_potential';
    fval_list=[fval_list;fval_potential];
    if ~isempty(con_list)
        con_list=[con_list;con_potential];con_potential=con_potential';
    end
    if ~isempty(coneq_list)
        coneq_list=[coneq_list;coneq_potential];coneq_potential=coneq_potential';
    end
    
    % updata penalty factor and lagrangian
    if expensive_nonlcon_flag
        lambda_con=lambda_con+2*miu*max(con_potential',-lambda_con/2/miu);
        lambda_coneq=lambda_coneq+2*miu*coneq_potential';
        if miu < miu_max
            miu=gama*miu;
        else
            miu=miu_max;
        end
    end
    
    % Step 5
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
        if (iteration > iteration_max*0.2 && ...
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
    
    % Step 6
    % TR guideline first iteration
    if iteration == 2
        bou_range_nomlz=0.5;
    else
        bou_range_nomlz_old=bou_range_nomlz;
        
        r=(fval_potential_old-fval_potential)/(fval_potential_old-fval_potential_predict);
        x_dis=norm((x_potential_old-x_potential)./(up_bou-low_bou),2);
        if x_dis <= range_min
            x_dis=0.1;
        end
        
        % scale only can be in (range_max-range_min)
        scale=enlarge_range*x_dis;
        % add rand influence avoid stable
        if scale > (range_max-range_min)
            scale=(range_max-range_min);
        end

        % mapping r into range_min - min(enlarge_range*x_dis,(range_max-range_min))
        % bou_range_nomlz=scale/2 while r=0
        bou_range_nomlz=scale/(1+exp(-(r-0.5)))+range_min;
        
        if abs(bou_range_nomlz-bou_range_nomlz_old) < torlance
           bou_range_nomlz=bou_range_nomlz*rand(); 
        end
    end
    bou_range=bou_range_nomlz.*(up_bou-low_bou);
    
    % updata trust range
    low_bou_iter=x_best-bou_range;
    low_bou_iter=max(low_bou_iter,low_bou);
    up_bou_iter=x_best+bou_range;
    up_bou_iter=min(up_bou_iter,up_bou);
    
    % Step 7
    % check whether exist data
    [x_list_exist,~,~,~]=dataLibraryLoad...
        (data_library_name,low_bou_iter,up_bou_iter);
    
    % generate latin hypercube sequence
    [~,x_updata_list,~]=getLatinHypercube...
        (sample_number,variable_number,x_list_exist,...
        low_bou_iter,up_bou_iter,cheapcon_function);
    
    x_potential_old=x_potential;
    fval_potential_old=fval_potential;
    fval_best_old=fval_best;
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
            (data_library_name,model_function,x_list,...
            x_exit_list,fval_exit_list,con_exit_list,coneq_exit_list,...
            low_bou,up_bou,protect_range)
        % function updata data with same_point_avoid protect
        % return fval
        % all list is x_number x variable_number matrix
        %
        NFE=0;
        x_updata_list=[];fval_updata_list=[];con_updata_list=[];coneq_updata_list=[];
        for x_index__=size(x_list,1)
            x_updata__=x_list(x_index__,:);
            
            % check x_potential if exist in data library
            % if not, updata data libraray
            distance__=sum(((x_updata__-x_exit_list)./(low_bou'-up_bou')).^2,2);
            [~,min_index__]=min(distance__);
            if distance__(min_index__) < protect_range^2
                x_updata__=x_exit_list(min_index__,:);
                fval_updata__=fval_exit_list(min_index__,:);
                con_updata__=[];
                coneq_updata__=[];
                if ~isempty(con_exit_list)
                    con_updata__=(con_exit_list(min_index__,:));
                end
                if ~isempty(coneq_exit_list)
                    coneq_updata__=(coneq_exit_list(min_index__,:));
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

function [x_best,fval_best]=findMinMSPFmincon...
    (object_function,x_initial,low_bou,up_bou,...
    cheapcon_function)
% find min fval use MSP guideline
% MSP: object_funtion is object_function (generate by surrogate model)
% nonlcon_function generate by surrogate model
% use ga as optimal method
%

% optiaml
fmincon_option=optimoptions('fmincon','display','none','algorithm','sqp');
[x_best,fval_best,~,~]=fmincon...
    (object_function,x_initial,[],[],[],[],low_bou,up_bou,cheapcon_function,fmincon_option);

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
max_cheapco_list=zeros(size(x_list,1),1);
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
        max_cheapco_list(x_index)=max_cheapco_list(x_index)+...
            sum(max(con,0))+sum(coneq.*coneq);
    end
end

con_judge_list=(max_nonlcon_list > nonlcon_torlance)+...
    (max_cheapco_list > 0);
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

% get initial X_new_nomalize by lhsdesign
X_new_nomlz=lhsdesign(x_new_number,variable_number);
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