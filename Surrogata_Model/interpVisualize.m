function interpVisualize...
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

axes_handle=axes(figure_handle);

x_list=model.X;
y_list=model.Y;
if iscell(x_list)
    x_list=x_list{1};
    y_list=y_list{1};
end
predict_function=model.predict_function;

% get boundary
if isempty(low_bou)
    low_bou=min(x_list,[],1);
end
if isempty(up_bou)
    up_bou=max(x_list,[],1);
end

if size(low_bou,1) ~= size(low_bou,1)
    error('interpolationRadialBasisVisualize: boundary incorrect');
end
if size(low_bou,1) > 2
    error('interpolationRadialBasisVisualize: dimension large than two');
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
    [X_draw,Y_draw]=meshgrid(low_bou(1):d_bou(1):(low_bou(1)+grid_number*d_bou(1)),...
        low_bou(2):d_bou(2):(low_bou(2)+grid_number*d_bou(2)));
    try
        % generate all predict list
        predict_x=[X_draw(:),Y_draw(:)];
        [predict_result]=predict_function(predict_x);
        predict_result=reshape(predict_result,grid_number+1,grid_number+1);
    catch
        % if not support multi input
        predict_result=zeros(grid_number+1);
        for x_index=1:grid_number+1
            for y_index=1:grid_number+1
                predict_x=([x_index,y_index]-1).*d_bou+low_bou;
                [predict_result(y_index,x_index)]=predict_function(predict_x);
            end
        end
    end
    surf(axes_handle,X_draw,Y_draw,predict_result,'FaceAlpha',0.5,'EdgeColor','none');
    line(axes_handle,x_list(:,1),x_list(:,2),y_list,'Marker','o','LineStyle','none');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(3);
end
end
