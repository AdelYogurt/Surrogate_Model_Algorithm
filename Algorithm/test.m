clc;
clear;
close all hidden;

load('matlab.mat');

tic;
for i = 1:10
[x_list,fval_list,con_list,coneq_list] = rankData...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance);
end
toc;

function [x_list,fval_list,con_list,coneq_list,vio_list] = rankData...
    (x_list,fval_list,con_list,coneq_list,...
    cheapcon_function,nonlcon_torlance)
% rank data base on feasibility rule
% infeasible is rank by sum of constraint
% torlance to cheapcon_function is 0
%
if nargin < 6 || isempty(nonlcon_torlance)
    nonlcon_torlance = 0;
end
if nargin < 5
    cheapcon_function = [];
end

[x_number,~] = size(x_list);
vio_list = zeros(x_number,1);
if ~isempty(con_list)
    vio_list = vio_list+sum(max(con_list-nonlcon_torlance,0),2);
end
if ~isempty(coneq_list)
    vio_list = vio_list+sum((abs(coneq_list)-nonlcon_torlance),2);
end

% add cheap con
for x_index = 1:size(x_list,1)
    if ~isempty(cheapcon_function)
        [con,coneq] = cheapcon_function(x_list(x_index,:));
        vio_list(x_index) = vio_list(x_index)+...
            sum(max(con,0))+sum(max(abs(coneq),0));
    end
end

% rank data
% infeasible data rank by violation, feasible data rank by fval
feasi_boolean_list = vio_list <= 0;
all = 1:x_number;
feasi_index_list = all(feasi_boolean_list);
infeasi_index_list = all(~feasi_boolean_list);
[~,index_list] = sort(fval_list(feasi_index_list));
feasi_index_list = feasi_index_list(index_list);
[~,index_list] = sort(vio_list(infeasi_index_list));
infeasi_index_list = infeasi_index_list(index_list);
index_list = [feasi_index_list,infeasi_index_list];

% rank by index_list
x_list = x_list(index_list,:);
fval_list = fval_list(index_list);
if ~isempty(con_list)
    con_list = con_list(index_list,:);
end
if ~isempty(coneq_list)
    coneq_list = coneq_list(index_list,:);
end
vio_list = vio_list(index_list);

end
