clc;
clear;
close all hidden;

sample_number=20;
variable_number=10;
low_bou=-3*ones(variable_number,1);
up_bou=3*ones(variable_number,1);
% x_exist_list=[0.1,0.5,0.3;0.6,0.3,0.6];
x_exist_list=[];

radius=2;

A=[];
B=[];
Aeq=[];
Beq=[];
cheapcon_function=@(x) sum(x.^2)-radius*radius;
% nomlz cheapcon_function
% optiaml
ga_option=optimoptions('ga','FunctionTolerance', 1e-2,...
    'PopulationSize',max(10,2*variable_number),...
    'MaxGenerations',50,...
    'display','none');
[~,fval_max]=ga...
    (@(x) -cheapcon_function(x),variable_number,A,B,Aeq,Beq,low_bou',up_bou',[],ga_option);
if fval_max < 0
    cheapcon_function=@(x) -cheapcon_function(x)/fval_max;
end

% cheapcon_function=[];

tic
[X,X_supply,distance_min_nomlz]=getLatinHypercube...
    (sample_number,variable_number,x_exist_list,...
    low_bou,up_bou)
toc
% scatter(X(:,1),X(:,2))
% hold on;
% rectangle('Position',[-radius,-radius,2*radius,2*radius],'Curvature',[1 1])
% hold off;

% bou=[low_bou,up_bou]';
% axis(bou(:));axis equal


function [X,X_new,distance_min_nomlz]=getLatinHypercube...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou)
% generate sample sequence latin hypercube
% iteration optimal method is used(Netwon method)
% not support cheap_function constraints
% sample number is total point in area
% default low_bou is 0, up_bou is 1, cheapcon_function is []
% low_bou and up_bou is colume vector
% x_exist_list, x_list, supply_x_list is row vector
% x_exist_list should meet bou
%
% Copyright 2022 Adel
%
if nargin < 5
    if nargin < 3
        if nargin < 2
            error('getLatinHypercube: lack variable_number');
        end
    end
    low_bou=zeros(variable_number,1);
    up_bou=ones(variable_number,1);
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
X_supply_quasi_nomlz=rand(x_initial_number,variable_number);

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

% x is nomalize, so constraint function should change
low_bou_nomlz=zeros(variable_number,1);
up_bou_nomlz=ones(variable_number,1);

iteration=1;
while iteration <= iteration_max
    for X_new_index=1:x_new_number
        % iteration
        object_function=@(x) objectFunctionXPlace...
            (x,[X_new_nomlz(1:X_new_index-1,:);X_new_nomlz(X_new_index+1:end,:);X_exist_nomlz],...
            low_bou_nomlz,up_bou_nomlz);
        %         drawFunction(object_function,low_bou_nomlz,up_bou_nomlz,100,-inf,1)
        %         fminunc_options=optimoptions('fminunc','display','none','StepTolerance',1e-2);
        %         [x,~,~,output]=fminunc(object_function,X_new(X_new_index,:)',fminunc_options);
        x=optimalNewton(object_function,X_new_nomlz(X_new_index,:)');
        x(find(x > 1))=1;
        x(find(x < 0))=0;
        X_new_nomlz(X_new_index,:)=x';
    end
    iteration=iteration+1;
end
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);
X_new=X_new_nomlz.*(up_bou'-low_bou')+low_bou';
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
                fval=zeros(grid_number+1,1);
                for x_index__=1:(grid_number+1)
                    fval(x_index__)=object_function(X__(x_index__));
                end
                line(axes_handle,X__,fval);
                xlabel('X');
                ylabel('value');
                
            case 2
                d_bou=(up_bou-low_bou)/grid_number;
                [X__,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
                fval=zeros(grid_number+1);
                for x_index__=1:grid_number+1
                    for y_index=1:grid_number+1
                        predict_x=([x_index__;y_index]-1).*d_bou+low_bou;
                        fval(x_index__,y_index)=object_function(predict_x);
                    end
                end
                fval(find(fval > Y_max))=Y_max;
                fval(find(fval < Y_min))=Y_min;
                surf(axes_handle,X__',Y',fval,'FaceAlpha',0.5,'EdgeColor','none');
                xlabel('X');
                ylabel('Y');
                zlabel('value');
        end
        
    end
    function x=optimalNewton(object_function,x,torlance,max_iteration__)
        % simple newton method optimal
        %
        if nargin < 4
            max_iteration__=50;
            if nargin < 3
                torlance=1e-3;
                if nargin < 2
                    error('getLatinHypercube:optimalNewton: lack x');
                end
            end
        end
        done__=0;
        iteration__=1;
        
        while ~done__
            [~,gradient,hessian]=object_function(x);
            %             if rcond(hessian) <1e-6
            %                 disp('error');
            %             end
            x=x-hessian\gradient;
            if iteration__ >= max_iteration__
                done__=1;
            end
            if norm(gradient,2) < torlance
                done__=1;
            end
            iteration__=iteration__+1;
        end
    end
    function [fval,gradient,hessian]=objectFunctionXPlace...
            (x,X_surplus,low_bou,up_bou)
        % function describe distance between X and X_supply
        % X__ is colume vector and X_supply__ is matrix which is num-1 x var
        % low_bou_limit__ and up_bou_limit__ is colume vector
        % variable in colume
        %
        [x_number,variable_number__]=size(X_surplus);
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
        
        function [gradient,hessian]=differ(differ_function,x,step,variable_number)
            % differ function to get gradient and hessian
            %
            if nargin < 4
                variable_number=size(x,1);
                if nargin < 3
                    step=1e-6;
                end
            end
            fval__=zeros(variable_number,3);
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
                fval__(variable_index__,3)=differ_function(x_forward__);
                
                gradient(variable_index__)=...
                    (fval__(variable_index__,3)-fval__(variable_index__,1))/2/step;
            end
            
            % hessian
            for variable_index__=1:variable_number
                hessian(variable_index__,variable_index__)=...
                    (fval__(variable_index__,3)-2*fval__(variable_index__,2)+fval__(variable_index__,1))/step/step;
                for variable_index_next__=variable_index__+1:variable_number
                    x_for_for=x;
                    x_for_for(variable_index__)=x_for_for(variable_index__)+step;
                    x_for_for(variable_index_next__)=x_for_for(variable_index_next__)+step;
                    fval_1__=differ_function(x_for_for);
                    
                    x_for_back=x;
                    x_for_back(variable_index__)=x_for_back(variable_index__)+step;
                    x_for_back(variable_index_next__)=x_for_back(variable_index_next__)-step;
                    fval_2__=differ_function(x_for_back);
                    
                    x_back_for=x;
                    x_back_for(variable_index__)=x_back_for(variable_index__)-step;
                    x_back_for(variable_index_next__)=x_back_for(variable_index_next__)+step;
                    fval_3__=differ_function(x_back_for);
                    
                    x_back_back=x;
                    x_back_back(variable_index__)=x_back_back(variable_index__)-step;
                    x_back_back(variable_index_next__)=x_back_back(variable_index_next__)-step;
                    fval_4__=differ_function(x_back_back);
                    
                    hessian(variable_index__,variable_index_next__)=(...
                        fval_1__-fval_2__-fval_3__+fval_4__...
                        )/4/step/step;
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
            % only search in min_distance
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