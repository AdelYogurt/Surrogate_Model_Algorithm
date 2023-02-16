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

radialbasis_model=interpolationRadialBasisPreModel(x_list,fval_list);
figure_handle=figure(1);
interpolationVisualize(radialbasis_model,low_bou,up_bou,figure_handle)

predict_function=radialbasis_model.predict_function;

function radialbasis_model=interpolationRadialBasisPreModel...
    (X,Y,basis_function)
% radial basis function interpolation pre model function version 1
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

[beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
    (X_nomlz,Y_nomlz,basis_function,x_number);

% initialization predict function
predict_function=@(predict_x) interpolationRadialBasisPredictor...
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

    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
            (X,Y,basis_function,x_number)
        % interpolation polynomial responed surface core function
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
    function [predict_y]=interpolationRadialBasisPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            beta,basis_function,predict_x)
        % radial basis function interpolation predict function
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
