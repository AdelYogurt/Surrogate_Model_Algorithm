clc;
clear;
close all hidden;

sample_number=12;
variable_number=2;
low_bou=-0*ones(variable_number,1);
up_bou=1*ones(variable_number,1);
% x_exist_list=[0.1,0.5;0.3,0.6];
x_exist_list=[];

cheapcon_function=@(x) sum(x.^2)-4;

tic;
[X,X_supply,distance_min_normalize]=getLatinHypercube...
    (sample_number,variable_number,x_exist_list,...
    low_bou,up_bou)
toc;

scatter(X(:,1),X(:,2))
% bou=[low_bou,up_bou]';
% axis(bou(:));

function [X,X_new,distance_min_normalize]=getLatinHypercube...
    (sample_number,variable_number,x_exist_list,...
    low_bou,up_bou,cheapcon_function)
% generate sample sequence latin hypercube
% DSLE method is used
% each dimension will be pick up with first dimension as space
% find max point in each space, construct as final result
% sample number is total point in area
% default low_bou is 0, up_bou is 1
% low_bou and up_bou is colume vector
% x in x_exist_list, x_list, supply_x_list is row vector
% x_exist_list should conform bou
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

% check x_exist_list if meet boundary
if ~isempty(x_exist_list)
    index=find(x_exist_list < low_bou');
    index=[index,find(x_exist_list > up_bou')];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list out of boundary');
    end
    if size(x_exist_list,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
end

x_new_number=sample_number-size(x_exist_list,1);
if x_new_number < 0
    X=x_exist_list;
    X_new=[];
    distance_min_normalize=getMinDistance((x_exist_list-low_bou')./(up_bou'-low_bou'));
    return;
end

% initialize
feasible_list=1:sample_number;
feasible_list=repmat(feasible_list,1,variable_number);

% process x exist, x exist correspond grid exist
if ~isempty(x_exist_list)
    x_exist_list_normalize=(x_exist_list-low_bou')./(up_bou'-low_bou');
    
    % change x_exist_list into grid
    grid_exist_list=zeros(size(x_exist_list_normalize,1),variable_number);
    for x_exist_index=1:size(x_exist_list_normalize,1)
        for variable_index=1:variable_number
            feasible=round(x_exist_list_normalize(x_exist_index,variable_index)*(sample_number-1))+1;
            
            % if feasible have exist in feasible_list, choose another grid
            if sum(grid_exist_list(:,variable_index)==feasible)
                practicable=1:sample_number;
                practicable(grid_exist_list(1:x_exist_index-1,variable_index))=[];
                distance=(practicable-feasible).^2;
                [~,index_min]=min(distance);
                feasible=practicable(index_min);
            end
            
            grid_exist_list(x_exist_index,variable_index)=feasible;
        end
    end
     
    % remove exist grid in feasiable range
    remove_index_list=grid_exist_list(:);
    place=repmat(0:variable_number-1,[size(x_exist_list_normalize,1),1]);
    place=place(:)*sample_number;
    feasible_list(remove_index_list+place)=[];
    
    if (length(feasible_list)/variable_number) ~= ...
            (sample_number-size(x_exist_list_normalize,1))
        save('matlab.mat');
        error('getLatinHypercube: x_exist_list dimension is repeat');
    end
end

% choose gird by max range
grid_list=zeros(sample_number,variable_number);
if ~isempty(x_exist_list)
   grid_list(1:size(x_exist_list,1),:)=grid_exist_list; 
end
for sample_new_index=size(x_exist_list,1)+1:sample_number
    grid_available_number=sample_number-sample_new_index+1;
    place_base=0:variable_number-1;
    place_base=place_base*grid_available_number;
    
    if sample_new_index==1
        % the first one gird can rand select
        grid_index_supply_list=randsrc(1,variable_number-1,1:sample_number);
        grid_index_list=[1,grid_index_supply_list];
    else
        % minimize the each grid to existing point distance
        % first colume is constraint min to max, donot select
        % because multiple variable number, grid_list available have to
        % be row list which length is grid_number^(variable_number-1)
        if variable_number==1
            grid_index_list=1;
            grid_list(sample_new_index,:)=feasible_list(grid_index_list);
        else
            grid_index_list=ones(1,variable_number);
            for variable_index=2:variable_number
                % find max distance in supply dimension
                grid_list__=[grid_list(:,1),grid_list(:,variable_index)];
                feasible_list__=[feasible_list(1:grid_available_number),...
                    feasible_list(((variable_index-1)*grid_available_number+1):variable_index*grid_available_number)];
                grid_index_list(variable_index)=...
                    getMinIndex(grid_list__,feasible_list__,...
                    grid_available_number,2,sample_new_index);
            end
        end
    end
    grid_list(sample_new_index,:)=feasible_list(grid_index_list+place_base);
    
    % remove index_list from feasiale list
    feasible_list(grid_index_list+place_base)=[];
end

% change gird into x
X_normalize=(grid_list-1)/(sample_number-1);
X_new_normalize=(grid_list(size(x_exist_list,1)+1:end,:)-1)/(sample_number-1);

distance_min_normalize=getMinDistance(X_normalize);
% normailze
X=X_normalize.*(up_bou'-low_bou')+low_bou';
X_new=X_new_normalize.*(up_bou'-low_bou')+low_bou';

    function grid_index_list=getMinIndex...
            (grid_list__,feasible_list__,...
        grid_available_number,variable_number,sample_new_index)
       % select grid for x_list, x_list is two dimension projection
       %
       distance_list=zeros(1,grid_available_number^(variable_number-1));
       
       for grid_index=1:length(distance_list)
           % get index_list from grid_index
           grid_index_list_supply__=translateIndex...
               (grid_index,grid_available_number,variable_number);
           
           % becasue feasible list is a row vector
           place_base_supply=1:variable_number-1;
           place_base_supply=place_base_supply*grid_available_number;
           
           % calculate grid to each exist point distance and find min distance
           grid_exist_number=sample_new_index-1;
           
           distance__=sum((grid_list__(1:grid_exist_number,2:variable_number)-...
               feasible_list__(grid_index_list_supply__+place_base_supply)).^2,2)+...
               (grid_list__(1:grid_exist_number,1)-...
               feasible_list__(1)).^2;
           distance_list(grid_index)=min(distance__);
       end
       [~,grid_index]=max(distance_list);
       grid_index=grid_index(1);
       grid_index_list=translateIndex...
           (grid_index,grid_available_number,variable_number);
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
    function index_list=translateIndex(index,grid_number,variable_number)
        % change index to list by grid_number and variable_number
        % index_list is 1 x (variable_number-1)
        %
        index_list=zeros(1,variable_number-1);
        for grid_list_index__=length(index_list):-1:1
            if grid_list_index__==1
                index_list(grid_list_index__)=index;
            else
                index_list(grid_list_index__)=ceil(index/...
                    grid_number^(grid_list_index__-1));
                index=index-...
                    grid_number^(grid_list_index__-1)*...
                    (index_list(grid_list_index__)-1);
            end
        end
    end
end