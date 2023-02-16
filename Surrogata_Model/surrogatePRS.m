clc;
clear;
close all hidden;

% flag=10;
% low_bou=[-3;-3];
% up_bou=[3;3];
% data=importdata('optimalSurrogateRespSurf_result.txt');
% x_list=data(1:flag,4:5);
% fval_list=data(1:flag,6);
% 
% respsurf_model=interpolationRespSurfPreModel(x_list,fval_list);
% 
% figure_handle=figure(1);
% interpolationRespSurfVisualize(respsurf_model,low_bou,up_bou,figure_handle)
% 
% X_new=data(flag+1:end,4:5);
% fval_new=data(flag+1:end,6);
% respsurf_model=interpolationRespSurfUpdata...
%     (respsurf_model,X_new,fval_new);
% 
% figure_handle=figure(2);
% interpolationRespSurfVisualize(respsurf_model,low_bou,up_bou,figure_handle)

% x_list=[2;3;4];
% fval_list=[2;3;4];
% low_bou=[2];
% up_bou=[4];

low_bou=[-3;-3];
up_bou=[3;3];
x_list=[1,2;
    1.5,-3;
    -2.2,1.5;
    -1,3;
    0.6,0.2;
    -0.8,0.1];
fval_list=zeros(6,1);
for x_index=1:size(fval_list,1)
    x=x_list(x_index,:);
    fval_list(x_index)=1+2*x(1)+3*x(2)+4*x(1)*x(1)+5*x(2)*x(2)+6*x(1)*x(2);
end

respsurf_model=interpolationRespSurfPreModel(x_list,fval_list);
interpolationRespSurfVisualize(respsurf_model,low_bou,up_bou)

function respsurf_model=interpolationRespSurfPreModel(X,Y)
% polynomial response surface interpolation pre model function
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
% beta is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
x_number=size(X,1);
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

beta=interpolationRespSurf(X_nomlz,Y_nomlz);

% initialization predict function
predict_function=@(predict_x) interpolationRespSurfPredictor...
    (aver_X,stdD_X,aver_Y,stdD_Y,beta,predict_x);

respsurf_model.X=X;
respsurf_model.Y=Y;
respsurf_model.X_normalize=X_nomlz;
respsurf_model.Y_normalize=Y_nomlz;
respsurf_model.beta=beta;
respsurf_model.aver_X=aver_X;
respsurf_model.stdD_X=stdD_X;
respsurf_model.aver_Y=aver_Y;
respsurf_model.stdD_Y=stdD_Y;
respsurf_model.predict_function=predict_function;

    function beta=interpolationRespSurf(X,Y)
        % interpolation polynomial responed surface core function
        % calculation beta
        %
        [x_number__,variable_number__]=size(X);
        
        X_inner__=zeros(x_number__,(variable_number__+1)*(variable_number__+2)/2);
        x_cross__=zeros(1,(variable_number__-1)*variable_number__/2);
        for x_index=1:x_number__
            x__=X(x_index,:);
            x_sqrt__=x__.^2;
            cross_index__=1;
            for i_index=1:variable_number__
                for j_index=i_index+1:variable_number__
                    x_cross__(cross_index__)=x__(i_index)*x__(j_index);
                    cross_index__=cross_index__+1;
                end
            end
            X_inner__(x_index,:)=[1,x__,x_sqrt__,x_cross__];
        end
        
        X_inter_X_inter__=X_inner__'*X_inner__;
        beta=X_inter_X_inter__\X_inner__'*Y;
    end
    function [predict_y]=interpolationRespSurfPredictor...
            (aver_X,stdD_X,aver_Y,stdD_Y,beta,predict_x)
        % polynomial response surface interpolation predict function
        % input predict_x and respsurf_model model
        % predict_x is row vector
        % output the predict value
        %
        % Copyright 2022 Adel
        %
        if size(predict_x,1) > 1
            predict_x=predict_x';
        end
        variable_number__=size(predict_x,2);
        
        % normalize data
        predict_x=(predict_x-aver_X)./stdD_X;
        
        % predict value
        x_cross__=zeros(1,(variable_number__-1)*variable_number__/2);
        x_sqrt__=predict_x.^2;
        cross_index__=1;
        for i_index=1:variable_number__
            for j_index=i_index+1:variable_number__
                x_cross__(cross_index__)=predict_x(i_index)*predict_x(j_index);
                cross_index__=cross_index__+1;
            end
        end
        x_inter=[1,predict_x,x_sqrt__,x_cross__];
        
        % predict variance
        predict_y=x_inter*beta;
        
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

axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end

x_list=model.X;
y_list=model.Y;
predict_function=model.predict_function;
[x_number,variable_number]=size(x_list);

% get boundary
if isempty(low_bou)
    low_bou=min(x_list,[],1)';
end
if isempty(up_bou)
    up_bou=max(x_list,[],1)';
end

% find point in boundary and label
index_list=[];
for x_index=1:x_number
    x=x_list(x_index,:)';
    if ~(sum(x > up_bou) || sum(x < low_bou))
        index_list=[index_list;x_index];
    end
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
    line(axes_handle,x_list(index_list,:),y_list(index_list,:),'Marker','o','LineStyle','none');
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
    line(axes_handle,x_list(index_list,1),x_list(index_list,2),y_list(index_list,:),'Marker','o','LineStyle','none');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
end
