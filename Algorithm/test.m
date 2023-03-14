clc;
clear;
close all hidden;

load('matlab.mat')

benchmark = BenchmarkFunction();

% variable_number = 20;
% object_function = @(x) benchmark.singleEP20Object(x);
% low_bou = ones(1,variable_number)*-5;
% up_bou = ones(1,variable_number)*5;

variable_number = 4;
nonlcon_function = @(x) benchmark.singlePVD4Nonlcon(x);
low_bou = [0,0,0,0];
up_bou = [1,1,50,240];

basis_function = @(r) r.^3;

% X = lhsdesign(50,variable_number).*(up_bou-low_bou)+low_bou;
% Y = nonlcon_function(X);

X = x_list;
% Y = con_nomlz_list(:,1);
Y = fval_nomlz_list;

[predict_function,radialbasis_model] = interpRadialBasisPreModel...
    (X,Y,basis_function);

% drawFunction(nonlcon_function,low_bou,up_bou,[],[],[],[],[3,4]);
% drawFunction(object_function,low_bou,up_bou,[],[],[],[],[1,2]);
drawFunction(predict_function,low_bou,up_bou,[],[],[],[],[1,2]);

% load('matlab.mat')
% 
% basis_function = @(r) r.^3;
% 
% interp_model_function = @(X,Y) interpRadialBasisPreModel(X,Y,basis_function);
% 
% x_number = size(x_list,1);
% 
% x_initial = x_best;
% distance = sum(((x_initial-x_list)./(up_bou-low_bou)).^2,2);
% [~,index_list] = sort(distance);
% 
% R2_list = zeros(x_number,1);
% for x_index = 5:x_number
%     RBF_number = x_index;
% 
%     index_local_list = index_list(1:RBF_number);
%     X = x_list(index_local_list,:);
%     Y = fval_list(index_local_list,:);
% 
% %     error_list = interpSampleVerify(interp_model_function,X,Y);
% 
%     [predict_function,radialbasis_model] = interpRadialBasisPreModel(X,Y,basis_function);
%     error_list = radialbasis_model.beta./diag(radialbasis_model.inv_radialbasis_matrix);
% 
%     R2_list(x_index) = 1-sum(error_list.^2)/sum((Y-mean(Y)).^2);
% end
% line(1:x_number,R2_list);

function error_list = interpSampleVerify(interp_model_function,X,Y)
% select each x in X as check point
% calculate predict error at each y
%
x_number = size(X,1);
error_list = zeros(x_number,1);
for x_index = 1:x_number
    [predict_function,~] = interp_model_function...
        (X([1:(x_index-1),(x_index+1):x_number],:),Y([1:x_index-1,x_index+1:x_number],:));
    error_list(x_index) = predict_function(X(x_index,:))-Y(x_index);
end
end

function [predict_function,radialbasis_model] = interpRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interp pre model function
% input initial data X,Y,which are real data
% X,Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model,include X,Y,base_function
% and predict_function
%
% Copyright 2023 Adel
%
if nargin < 3
    basis_function = [];
end

[x_number,variable_number] = size(X);

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__ = find(stdD_X == 0);
if ~isempty(index__),stdD_X(index__) = 1;end
index__ = find(stdD_Y == 0);
if ~isempty(index__),stdD_Y(index__) = 1;end
X_nomlz = (X-aver_X)./stdD_X;
Y_nomlz = (Y-aver_Y)./stdD_Y;

if isempty(basis_function)
    c = (prod(max(X_nomlz)-min(X_nomlz))/x_number)^(1/variable_number);
    basis_function = @(r) exp(-(r.^2)/c);
end

% initialization distance of all X
X_dis = zeros(x_number,x_number);
for variable_index = 1:variable_number
    X_dis = X_dis+(X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end
X_dis = sqrt(X_dis);

[beta,rdibas_matrix,inv_rdibas_matrix] = interpRadialBasis...
    (X_dis,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function = @(X_predict) interpRadialBasisPredictor...
    (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,beta,basis_function);

radialbasis_model.X = X;
radialbasis_model.Y = Y;
radialbasis_model.radialbasis_matrix = rdibas_matrix;
radialbasis_model.inv_radialbasis_matrix=inv_rdibas_matrix;
radialbasis_model.beta = beta;

radialbasis_model.aver_X = aver_X;
radialbasis_model.stdD_X = stdD_X;
radialbasis_model.aver_Y = aver_Y;
radialbasis_model.stdD_Y = stdD_Y;
radialbasis_model.basis_function = basis_function;

radialbasis_model.predict_function = predict_function;

% abbreviation:
% num: number,pred: predict,vari: variable
    function [beta,rdibas_matrix,inv_rdibas_matrix] = interpRadialBasis...
            (X_dis,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);
        
        % stabilize matrix
        rdibas_matrix = rdibas_matrix+eye(x_number)*1e-6;

        % get inverse matrix
        inv_rdibas_matrix = rdibas_matrix\eye(x_number);
        
        % solve beta
        beta = inv_rdibas_matrix*Y;
    end

    function [Y_pred] = interpRadialBasisPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,beta,basis_function)
        % radial basis function interpolation predict function
        %
        [x_pred_num,~] = size(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred-aver_X)./stdD_X;
        
        % calculate distance
        X_dis_pred = zeros(x_pred_num,x_num);
        for vari_index = 1:vari_num
            X_dis_pred = X_dis_pred+...
                (X_pred_nomlz(:,vari_index)-X_nomlz(:,vari_index)').^2;
        end
        X_dis_pred = sqrt(X_dis_pred);
        
        % predict variance
        Y_pred = basis_function(X_dis_pred)*beta;
        
        % normalize data
        Y_pred = Y_pred*stdD_Y+aver_Y;
    end

end
