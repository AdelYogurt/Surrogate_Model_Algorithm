clc;
clear
close all hidden;

load('ZDT1.mat');

[~,index_list]=sort(fval_pareto(:,1));
fval_pareto=fval_pareto(index_list,:);
x_pareto=x_pareto(index_list,:);

line(fval_pareto(:,1),fval_pareto(:,2),'Marker','o','LineStyle','none');

real_pareto=[(0:0.01:1)',(1-(0:0.01:1).^0.5)'];
% real_pareto=[(0:0.01:1)',(1-(0:0.01:1).^2)'];
% real_pareto=[(0:0.01:1)',(1-(0:0.01:1).^0.5-(0:0.01:1).*sin(10*pi*(0:0.01:1)))'];pareto_index_list=getParetoFront(real_pareto,[],[]);real_pareto=real_pareto(pareto_index_list,:);

line(real_pareto(:,1),real_pareto(:,2));

IGD = mean(min(pdist2(fval_pareto,real_pareto),[],2))

function [pareto_index_list,feasible_index_list]=getParetoFront...
    (fval_nomlz_list,con_nomlz_list,coneq_nomlz_list,...
    pareto_torlance)
% distinguish pareto front of exist point
% use initial pareto definition(contrain compare first and fval compare)
%
if nargin < 4
    pareto_torlance=0;
end

pareto_index_list=[];% index of filter point list
feasible_index_list=[];% feasible point list

% select no domain filter
for fval_index=1:size(fval_nomlz_list,1)
    con_nomlz=[];
    if ~isempty(con_nomlz_list)
        con_nomlz=max(con_nomlz_list(fval_index,:));
    end
    coneq_nomlz=[];
    if ~isempty(coneq_nomlz_list)
        coneq_nomlz=max(abs(coneq_nomlz_list(fval_index,:)));
    end
    con_nomlz_max=max([con_nomlz;coneq_nomlz;0]); % con is large than 0 or 0
    
    if (con_nomlz_max <= pareto_torlance)
        feasible_index_list=[feasible_index_list;fval_index];
    end
    
    % notice con_nomlz_max is greater than or equal to zero
    % notice con_nomlz_pareto_max is greater than or equal to zero
    pareto_index=1;
    add_filter_flag=1;
    while pareto_index <= length(pareto_index_list)
        % compare x with exit pareto front point
        x_pareto_index=pareto_index_list(pareto_index,:);
        
        % contain constraint of x_filter
        con_pareto_nomlz=[];
        if ~isempty(con_nomlz_list)
            con_pareto_nomlz=max(con_nomlz_list(x_pareto_index,:));
        end
        coneq_pareto_nomlz=[];
        if ~isempty(coneq_nomlz_list)
            coneq_pareto_nomlz=max(coneq_nomlz_list(x_pareto_index,:));
        end
        con_nomlz_pareto_max=max([con_pareto_nomlz;coneq_pareto_nomlz;0]);
        
        % compare x with x_pareto
        % if this x is domain by pareto point, reject it
        if con_nomlz_max > con_nomlz_pareto_max
            % x constrain is large than pareto
            add_filter_flag=0;
            break;
        elseif con_nomlz_max <= pareto_torlance && con_nomlz_pareto_max <= pareto_torlance
            % both is feasiable point
            judge=fval_nomlz_list(fval_index,:)+pareto_torlance >= fval_nomlz_list(x_pareto_index,:);
            if ~sum(~judge)
                add_filter_flag=0;
                break;
            end
        end
        
        % if better or equal than exit pareto point, reject pareto point
        if con_nomlz_pareto_max > con_nomlz_max
            % pareto constrain is large than x
            pareto_index_list(pareto_index)=[];
            pareto_index=pareto_index-1;
        elseif con_nomlz_max <= pareto_torlance && con_nomlz_pareto_max <= pareto_torlance
            judge=fval_nomlz_list(fval_index,:) <= fval_nomlz_list(x_pareto_index,:)+pareto_torlance;
            if ~sum(~judge)
                pareto_index_list(pareto_index)=[];
                pareto_index=pareto_index-1;
            end
        end
        
        pareto_index=pareto_index+1;
    end
    
    % add into pareto list if possible
    if add_filter_flag
        pareto_index_list=[pareto_index_list;fval_index];
    end
end
end