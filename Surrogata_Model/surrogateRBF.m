clc;
clear;
close all hidden;


% data=importdata('optimalSurrogateRadialBasis_result.txt');
% x_list=data(1:6,end-1:end);
% fval_list=data(1:6,1);
% low_bou=[-3;-3];
% up_bou=[3;3];

load('matlab')

% x_list=[2;3;4];
% fval_list=[2;3;4];
% low_bou=[2];
% up_bou=[4];

radialbasis_model=interpRadialBasisPreModel(x_list,fval_list);
figure_handle=figure(1);
interpVisualize(radialbasis_model,low_bou,up_bou,figure_handle)

predict_function=radialbasis_model.predict_function;

function radialbasis_model=interpRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interp pre model function version 1
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
% beta is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
if nargin < 9
    X_nomlz=[];
    if nargin < 3
        basis_function=[];
    end
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
X_nomlz = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_nomlz = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

if isempty(basis_function)
    c=(prod(max(X_nomlz)-min(Y_nomlz))/x_number)^(1/variable_number);
    basis_function=@(r) exp(-(r'*r)/c);
%     basis_function=@(r) sqrt(r'*r+c*c);
end

[beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
    (X_nomlz,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function=@(predict_x) interpRadialBasisPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    beta,basis_function,predict_x);

radialbasis_model.X=X;
radialbasis_model.Y=Y;
radialbasis_model.X_normalize=X_nomlz;
radialbasis_model.Y_normalize=Y_nomlz;
radialbasis_model.radialbasis_matrix=rdibas_matrix;
radialbasis_model.inv_radialbasis_matrix=inv_rdibas_matrix;
radialbasis_model.beta=beta;

radialbasis_model.aver_X=aver_X;
radialbasis_model.stdD_X=stdD_X;
radialbasis_model.aver_Y=aver_Y;
radialbasis_model.stdD_Y=stdD_Y;
radialbasis_model.basis_function=basis_function;

radialbasis_model.predict_function=predict_function;

    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
            (X,Y,basis_function,x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix=zeros(x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                rdibas_matrix(rank_index__,colume_index__)=...
                    rdibas_matrix(colume_index__,rank_index__);
            end
            for colume_index__=rank_index__:x_number
                rdibas_matrix(rank_index__,colume_index__)=...
                    basis_function(X(rank_index__,:)'-X(colume_index__,:)');
            end
        end
        
        % stabilize matrix
        rdibas_matrix=rdibas_matrix+eye(x_number)*1e-6;
        
        inv_rdibas_matrix=inv(rdibas_matrix);
        beta=inv_rdibas_matrix*Y;
    end
    function [predict_y]=interpRadialBasisPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            beta,basis_function,predict_x)
        % radial basis function interp predict function
        % input predict_x and radialbasis_model model
        % predict_x is row vector
        % output the predict value
        %
        % Copyright 2022 Adel
        %
        predict_x=predict_x(:);
        
        [x_number__,~]=size(X_nomlz);
        
        % normalize data
        predict_x=(predict_x-aver_X')./stdD_X';
        
        % predict value
        X_inner_product__=zeros(x_number__,1);
        for index_i=1:x_number__
            X_inner_product__(index_i,1)=...
                basis_function(X_nomlz(index_i,:)'-predict_x);
        end
        
        % predict variance
        predict_y=beta'*X_inner_product__;
        
        % normalize data
        predict_y=predict_y*stdD_Y+aver_Y;
    end
end
