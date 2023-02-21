clc;
clear;
close all hidden;

sample_number=12;
variable_number=2;
low_bou=-0*ones(variable_number,1);
up_bou=1*ones(variable_number,1);

cheapcon_function=@(x) sum(x.^2)-4;

tic;
[X,X_supply,distance_min_nomlz]=getLatinHypercube...
    (sample_number,variable_number,[],...
    low_bou,up_bou,[])
toc;

scatter(X(:,1),X(:,2))
bou=[low_bou,up_bou]';
axis(bou(:));

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
    index=find(X_exist < low_bou');
    index=[index,find(X_exist > up_bou')];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list range error');
    end
    X_exist_nomlz=(X_exist-low_bou')./(up_bou'-low_bou');
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
            x_supply=X_supply_initial_nomlz(x_index,:).*(up_bou'-low_bou')+low_bou';
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
X_new=X_new_nomlz.*(up_bou'-low_bou')+low_bou';
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