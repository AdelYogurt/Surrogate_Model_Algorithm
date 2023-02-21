clc;
clear;
close all hidden;

sample_number=12;
variable_number=2;
low_bou=-0*ones(variable_number,1);
up_bou=1*ones(variable_number,1);
% x_exist_list=[0.5,0.5;0.3,0.6];
x_exist_list=[];

cheapcon_function=@(x) sum(x.^2)-4;

tic;
[X,X_supply,distance_min_normalize]=getLatinHypercube...
    (sample_number,variable_number,x_exist_list,...
    low_bou,up_bou,cheapcon_function)
toc;

scatter(X(:,1),X(:,2))
bou=[low_bou,up_bou]';
axis(bou(:));

function [X,X_new,distance_min_normalize]=getLatinHypercube...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou,cheapcon_function)
% generate sample sequence latin hypercube
% successive local enumeration method is used
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
if ~isempty(X_exist)
    index=find(X_exist < low_bou');
    index=[index,find(X_exist > up_bou')];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list out of boundary');
    end
    if size(X_exist,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
    X_exist_nomlz=(X_exist-low_bou')./(up_bou'-low_bou');
else
    X_exist_nomlz=[];
end

x_new_number=sample_number-size(X_exist,1);
if x_new_number < 0
    X=X_exist;
    X_new=[];
    distance_min_normalize=getMinDistance(X_exist_nomlz);
    return;
end

% initialize, feasible_list is usable grid list
% feasible_list is 1 x usable_grid_number*variable_number matrix
feasible_list=1:sample_number;
feasible_list=repmat(feasible_list,1,variable_number);

% process x exist, x exist correspond grid exist
if ~isempty(X_exist)
    X_exist_nomlz=(X_exist-low_bou')./(up_bou'-low_bou');
    
    % change x_exist_list into grid
    grid_exist_list=zeros(size(X_exist_nomlz,1),variable_number);
    for x_exist_index=1:size(X_exist_nomlz,1)
        for variable_index=1:variable_number
            % mapping x_exist to grid
            feasible=round(X_exist_nomlz(x_exist_index,variable_index)*(sample_number-1))+1;
            
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
    place=repmat(0:variable_number-1,[size(X_exist_nomlz,1),1]);
    place=place(:)*sample_number;
    feasible_list(remove_index_list+place)=[];
    
    if (length(feasible_list)/variable_number) ~= ...
            (sample_number-size(X_exist_nomlz,1))
        save('matlab.mat');
        error('getLatinHypercube: x_exist_list dimension is repeat');
    end
end

% choose gird by max range
grid_list=zeros(sample_number,variable_number);
if ~isempty(X_exist)
   grid_list(1:size(X_exist,1),:)=grid_exist_list; 
end
for sample_new_index=size(X_exist,1)+1:sample_number
    grid_available_number=sample_number-sample_new_index+1;
    place_base=0:variable_number-1;
    place_base=place_base*grid_available_number;
    
    if sample_new_index==1
        % the first one gird can rand select
        index_list=[1,randsrc(1,variable_number-1,1:sample_number)];
    else
        % first colume is constraint min to max, donot select
        % because multiple variable number, grid_list available have to
        % be row list which length is grid_number^(variable_number-1)
        % minimize the each grid to existing point distance
        % first find out each grid's min distance to exist point
        % then select min distance grid form list
        % index_list is the select grid index in feasible_list
        if variable_number==1
            index_list=1;
            grid_list(sample_new_index,:)=feasible_list(index_list);
        else
            % get all grid list
            all_grid=getAllGrid...
                (feasible_list,grid_available_number,variable_number);
            
            % calculate all grid list to exist point distance
            grid_exit_dis=zeros(grid_available_number^(variable_number-1),sample_new_index-1);
            for x_exist_index=1:sample_new_index-1
                grid_exit_dis(:,x_exist_index)=sum((grid_list(x_exist_index,:)-all_grid).^2,2);
            end
            
            [~,grid_index]=max(min(grid_exit_dis,[],2));
            grid_index=grid_index(1);
            grid_list(sample_new_index,:)=all_grid(grid_index,:);
            
            % transform grid index to index_list
            index_list=ones(1,variable_number);
            grid_index=grid_index-1;
            bit_index=variable_number;
            while(grid_index > 0)
                index_list(bit_index)=mod(grid_index,grid_available_number)+1;
                grid_index=floor(grid_index/grid_available_number);
                bit_index=bit_index-1;
            end
            
        end
    end
    
    % remove index_list from feasiale list
    feasible_list(index_list+place_base)=[];
end

% change gird into x
X_normalize=(grid_list-1)/(sample_number-1);
X_new_normalize=(grid_list(size(X_exist,1)+1:end,:)-1)/(sample_number-1);

distance_min_normalize=getMinDistance(X_normalize);
% normailze
X=X_normalize.*(up_bou'-low_bou')+low_bou';
X_new=X_new_normalize.*(up_bou'-low_bou')+low_bou';

    function all_grid=getAllGrid...
            (feasible_list,grid_available_number,variable_number)
        % get all grid list
        %
        all_grid=zeros(grid_available_number^(variable_number-1),variable_number);
        for variable_index__=2:variable_number
            g=ones(grid_available_number^(variable_number-variable_index__),1)*...
                feasible_list(1,((variable_index__-1)*grid_available_number+1):(variable_index__*grid_available_number));
            g=g(:);
            t=g*ones(1,grid_available_number^(variable_index__-2));
            all_grid(:,variable_index__)=t(:);
        end
        all_grid(:,1)=feasible_list(1);
    end
    function distance_min__=getMinDistance(x_list__)
        % get distance min from x_list
        % all x will be calculate
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
end