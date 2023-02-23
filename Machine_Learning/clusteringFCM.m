clc;
clear;
close all hidden;

x_number=50;
variable_number=2;

classify_number=2;

low_bou=ones(variable_number,1)*0;
up_bou=ones(variable_number,1);
m=2;

X=rand(x_number,variable_number);
FC_model=classifyFuzzyClustering(X,classify_number,low_bou,up_bou,m);

X=FC_model.X;
center_list=FC_model.center_list;

scatter(X(:,1),X(:,2));
line(center_list(:,1),center_list(:,2),'Marker','o','LineStyle','None','Color','r');

function FC_model=classifyFuzzyClustering...
    (X,classify_number,low_bou,up_bou,m)
% get fuzzy cluster model
% X is x_number x variable_number matrix
% center_list is classify_number x variable_number matrix
%
if nargin < 4
    up_bou=[];
    if nargin < 3
        low_bou=[];
    end
end
iteration_max=100;
torlance=1e-3;

[x_number,variable_number]=size(X);

% nomalization data
if isempty(low_bou)
    low_bou=min(X,[],1)';
end
if isempty(up_bou)
    up_bou=max(X,[],1)';
end
X_nomlz=(X-low_bou')./(up_bou'-low_bou');

% if x_number equal 1, clustering cannot done
if x_number==1
    FC_model.X=X;
    FC_model.X_normalize=X_nomlz;
    FC_model.center_list=X;
    FC_model.fval_loss_list=[];
    return;
end

U=zeros(classify_number,x_number);
center_list=rand(classify_number,variable_number)*0.5;
iteration=0;
done=0;
fval_loss_list=zeros(iteration_max,1);

% get X_center_dis_sq
X_center_dis_sq=zeros(classify_number,x_number);
for classify_index=1:classify_number
    for x_index=1:x_number
        X_center_dis_sq(classify_index,x_index)=...
            getSq((X_nomlz(x_index,:)-center_list(classify_index,:)));
    end
end

while ~done
    % updata classify matrix U
    for classify_index=1:classify_number
        for x_index=1:x_number
            U(classify_index,x_index)=...
                1/sum((X_center_dis_sq(classify_index,x_index)./X_center_dis_sq(:,x_index)).^(1/(m-1)));
        end
    end
    
    % updata center_list
    center_list_old=center_list;
    for classify_index=1:classify_number
        center_list(classify_index,:)=...
            sum((U(classify_index,:)').^m.*X_nomlz,1)./...
            sum((U(classify_index,:)').^m,1);
    end
    
    % updata X_center_dis_sq
    X_center_dis_sq=zeros(classify_number,x_number);
    for classify_index=1:classify_number
        for x_index=1:x_number
            X_center_dis_sq(classify_index,x_index)=...
                getSq((X_nomlz(x_index,:)-center_list(classify_index,:)));
        end
    end
    
%     plot(center_list(:,1),center_list(:,2));
    
    % forced interrupt
    if iteration > iteration_max
        done=1;
    end
    
    % convergence judgment
    if sum(sum(center_list_old-center_list).^2)<torlance
        done=1;
    end
    
    iteration=iteration+1;
    fval_loss_list(iteration)=sum(sum(U.^m.*X_center_dis_sq));
end
fval_loss_list(iteration+1:end)=[];
center_list=center_list.*(up_bou'-low_bou')+low_bou';

FC_model.X=X;
FC_model.X_normalize=X_nomlz;
FC_model.center_list=center_list;
FC_model.fval_loss_list=fval_loss_list;

    function sq=getSq(dx)
        % dx is 1 x variable_number matrix
        %
        sq=dx*dx';
    end
end
