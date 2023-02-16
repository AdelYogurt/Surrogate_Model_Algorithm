clc;
clear;
close all hidden;

% data=importdata('optimalSurrogateRadialBasis_result.txt');
% X=data(1:6,end-1:end);
% Y=data(1:6,1);
% low_bou=[-3;-3];
% up_bou=[3;3];

% load('matlab')

% X=[0.2;0.3;0.4];
% Y=[0.2;0.3;0.4];
% low_bou=0.2;
% up_bou=0.4;

low_bou=[-3;-3];
up_bou=[3;3];
x_list=getLatinHypercube(10,2,[],low_bou,up_bou);
fval_list=zeros(size(x_list,1),1);
for x_index=1:size(x_list,1)
    x=x_list(x_index,:);
    fval_list(x_index)=functionPKObject(x);
end

tic
ensemleradialbasis_model=surrogateEnsemleRadialBasisPreModel...
    (x_list,fval_list);
interpolationVisualize(ensemleradialbasis_model,low_bou,up_bou)
toc

function ensemleradialbasis_model=surrogateEnsemleRadialBasisPreModel...
    (X,Y)
% get ensemle radial basis function interpolation model function version 1
% use cubic interpolation optimal to decrese time use
% input initial data X, Y, which are real data
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% output is a radial basis model, include X, Y, base_function
% and predict_function
% beta is normalizede, so predict y is normalizede
%
% Copyright 2022 Adel
%
[x_number,variable_number]=size(X);

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index = find(stdD_X == 0);
if ~isempty(index),stdD_X(index)=1;end
index = find(stdD_Y == 0);
if ~isempty(index),stdD_Y(index)=1;end
X_nomlz = (X - repmat(aver_X,x_number,1)) ./ repmat(stdD_X,x_number,1);
Y_nomlz = (Y - repmat(aver_Y,x_number,1)) ./ repmat(stdD_Y,x_number,1);

c_initial=(prod(max(X_nomlz)-min(Y_nomlz))/x_number)^(1/variable_number);

% initial sq of x and x
x_sq_matrix=zeros(x_number);
for rank_index=1:x_number
    for colume_index=1:rank_index-1
        x_sq_matrix(rank_index,colume_index)=...
            x_sq_matrix(colume_index,rank_index);
    end
    for colume_index=rank_index:x_number
        x_sq_matrix(rank_index,colume_index)=...
            sum((X_nomlz(rank_index,:)'-X_nomlz(colume_index,:)').^2);
    end
end

option_fmincon=optimoptions('fmincon','display','none','TolFun',1e-6,...
    'SpecifyObjectiveGradient',true);

% linear kernal function
basis_function_linear=@(x_sq,c) sqrt(x_sq)+c;
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) ones(x_number);
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_linear,c,rdibas_matrix_gradient_function);
[linear_c_backward,fval_backward,NFE]=optimalCubicInterpolation...
    (object_function,-1e2,-1e2,1e2,1e-3);
[linear_c_forward,fval_forward,NFE]=optimalCubicInterpolation...
    (object_function,1e2,-1e2,1e2,1e-3);
if fval_forward < fval_backward
    c_linear=linear_c_forward;
else
    c_linear=linear_c_backward;
end

% gauss kernal function
basis_function_gauss=@(x_sq,c) exp(-c*x_sq);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) -x_sq_matrix.*rdibas_matrix;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_gauss,c,rdibas_matrix_gradient_function);
[c_gauss,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% spline kernal function
basis_function_spline=@(x_sq,c) x_sq*log(x_sq*c+1e-3);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) x_sq_matrix.*x_sq_matrix./(x_sq_matrix*c+1e-3);
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_spline,c,rdibas_matrix_gradient_function);
[c_spline,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% triple kernal function
basis_function_triple=@(x_sq,c) (sqrt(x_sq)+c)^3;
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) 3*(sqrt(x_sq_matrix)+c).^3;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_triple,c,rdibas_matrix_gradient_function);
[c_triple,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% multiquadric kernal function
basis_function_multiquadric=@(x_sq,c) sqrt(x_sq+c);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) 0.5./rdibas_matrix;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_multiquadric,c,rdibas_matrix_gradient_function);
[c_binomial,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% inverse multiquadric kernal function
basis_function_inverse_multiquadric=@(x_sq,c) 1/sqrt(x_sq+c);
rdibas_matrix_gradient_function=@(x_number,x_sq_matrix,rdibas_matrix,c) -0.5*rdibas_matrix.^3;
object_function=@(c) objectFunctionRadiabasis....
    (x_sq_matrix,Y_nomlz,x_number,basis_function_inverse_multiquadric,c,rdibas_matrix_gradient_function);
% c_initial=1;
% [fval,gradient]=object_function(c_initial)
% [fval_diff,gradient_diff]=differ(object_function,c_initial)
% drawFunction(object_function,1e-1,10);
[c_inverse_binomial,~,NFE,~]=optimalCubicInterpolation...
    (object_function,c_initial,1e-2,1e2,1e-3);

% generate total model
basis_function_list={basis_function_linear;
    basis_function_gauss;
    basis_function_spline;
    basis_function_triple;
    basis_function_multiquadric;
    basis_function_inverse_multiquadric;};
c_list=[c_linear;c_gauss;c_spline;c_triple;c_binomial;c_inverse_binomial];

model_number=size(basis_function_list,1);
beta_list=zeros(x_number,1,model_number);
rdibas_matrix_list=zeros(x_number,x_number,model_number);
inv_rdibas_matrix_list=zeros(x_number,x_number,model_number);

% calculate model matrix and error
model_error_list=zeros(model_number,x_number);
for model_index=1:model_number
    basis_function=basis_function_list{model_index};
    c=c_list(model_index);
    [beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
        (x_sq_matrix,Y_nomlz,x_number,basis_function,c);
    beta_list(:,:,model_index)=beta;
    rdibas_matrix_list(:,:,model_index)=rdibas_matrix;
    inv_rdibas_matrix_list(:,:,model_index)=inv_rdibas_matrix;
    
    model_error_list(model_index,:)=(beta./...
        diag(inv_rdibas_matrix))';
end

% calculate weight of each model
C=model_error_list*model_error_list';
eta=trace(C)/x_number;
w=(C+eta*eye(model_number))\ones(model_number,1)/...
    (ones(1,model_number)/(C+eta*eye(model_number))*ones(model_number,1));
while min(w) < -0.05
    eta=eta*10;
    w=(C+eta*eye(model_number))\ones(model_number,1)/...
        (ones(1,model_number)/(C+eta*eye(model_number))*ones(model_number,1));
end

% initialization predict function
predict_function=@(predict_x) interpolationEnsemleRadialBasisPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    model_number,beta_list,basis_function_list,c_list,...
    w,predict_x);

ensemleradialbasis_model.X=X;
ensemleradialbasis_model.Y=Y;
ensemleradialbasis_model.X_normalize=X_nomlz;
ensemleradialbasis_model.Y_normalize=Y_nomlz;
ensemleradialbasis_model.aver_X=aver_X;
ensemleradialbasis_model.stdD_X=stdD_X;
ensemleradialbasis_model.aver_Y=aver_Y;
ensemleradialbasis_model.stdD_Y=stdD_Y;

ensemleradialbasis_model.basis_function_list=basis_function_list;
ensemleradialbasis_model.c_list=c_list;
ensemleradialbasis_model.beta_list=beta_list;
ensemleradialbasis_model.rdibas_matrix_list=rdibas_matrix_list;
ensemleradialbasis_model.inv_rdibas_matrix_list=inv_rdibas_matrix_list;
ensemleradialbasis_model.model_error_list=model_error_list;
ensemleradialbasis_model.w=w;

ensemleradialbasis_model.predict_function=predict_function;

    function [fval,gradient]=objectFunctionRadiabasis....
            (x_sq_matrix,Y,x_number,basis_function,c,rdibas_matrix_gradient_function)
        % MSE_CV function, simple approximation to RMS
        % basis_function input is c and x_sq
        %
        [beta__,rdibas_matrix__,inv_rdibas_matrix__]=interpolationRadialBasis...
            (x_sq_matrix,Y,x_number,basis_function,c);
        fval=0;
        U=zeros(x_number,1);
        for x_index=1:x_number
            U(x_index)=...
                beta__(x_index)/inv_rdibas_matrix__(x_index,x_index);
            fval=fval+(U(x_index))^2;
        end
        
        % calculate gradient
        if nargout > 1
            inv_rdibas_matrix_gradient=-inv_rdibas_matrix__*...
                rdibas_matrix_gradient_function...
                (x_number,x_sq_matrix,rdibas_matrix__,c)*inv_rdibas_matrix__;
            U_gradient=zeros(x_number,1);
            I=eye(x_number);
            for x_index=1:x_number
                U_gradient(x_index)=(I(x_index,:)*inv_rdibas_matrix_gradient*Y)/...
                    inv_rdibas_matrix__(x_index,x_index)-...
                    beta__(x_index)*(I(x_index,:)*inv_rdibas_matrix_gradient*I(:,x_index))/...
                    inv_rdibas_matrix__(x_index,x_index)^2;
            end
            
            gradient=2*sum(U.*U_gradient);
        end
    end
    function [beta,rdibas_matrix,inv_rdibas_matrix]=interpolationRadialBasis...
            (x_sq_matrix,Y,x_number,basis_function,c)
        % interpolation polynomial responed surface core function
        % calculation beta
        rdibas_matrix=zeros(x_number);
        for rank_index__=1:x_number
            for colume_index__=1:rank_index__-1
                rdibas_matrix(rank_index__,colume_index__)=...
                    rdibas_matrix(colume_index__,rank_index__);
            end
            for colume_index__=rank_index__:x_number
                rdibas_matrix(rank_index__,colume_index__)=...
                    basis_function(x_sq_matrix(rank_index__,colume_index__),c);
            end
        end
        
        % stabilize matrix
        rdibas_matrix=rdibas_matrix+eye(x_number)*1e-6;
        
        if rcond(rdibas_matrix) < eps || isnan(rdibas_matrix(1))
            disp('error');
        end
        
        inv_rdibas_matrix=inv(rdibas_matrix);
        beta=inv_rdibas_matrix*Y;
    end
    function [x_best,favl_best,NFE,output]=optimalCubicInterpolation...
            (object_function,x_initial,low_bou,up_bou,torlance,iteration_max)
        % cubic interpolation optimization, should provide fval and gradient
        % only work for one best(convex)
        %
        if nargin < 6
            iteration_max=[];
            if nargin < 5
                torlance=[];
                if nargin < 4
                    up_bou=[];
                    if nargin < 3
                        low_bou=[];
                        if nargin < 2
                            error('lack x initial');
                        end
                    end
                end
            end
        end
        
        INFORMATION_FLAG=0; % whether show information
        
        draw_range=0.001;
        draw_interval=draw_range*0.02;
        DRAW_FLAG=0;
        
        if isempty(iteration_max)
            iteration_max=10*length(x_initial);
        end
        
        if isempty(torlance)
            torlance=1e-6;
        end
        
        x=x_initial;
        done=0;
        iteration=0;
        NFE=0;
        result_x_list=[];
        result_fval_list=[];
        
        % decide which turn to search
        [fval,gradient]=object_function(x);NFE=NFE+1;
        result_x_list=[result_x_list;x];
        result_fval_list=[result_fval_list;fval];
        if gradient < -torlance
            direction=1;
        elseif gradient > torlance
            direction=-1;
        else
            done=1;
            x_best=x;
            favl_best=fval;
        end
        
        x_old=x;
        fval_old=fval;
        gradient_old=gradient;
        iteration=iteration+1;
        
        % move forward to first point
        if ~done
            x=x_old+direction*0.01;
            if x > up_bou
                x=up_bou;
            elseif x < low_bou
                x=low_bou;
            end
            [fval,gradient]=object_function(x);NFE=NFE+1;
            result_x_list=[result_x_list;x];
            result_fval_list=[result_fval_list;fval];
            quit_flag=judgeQuit...
                (x,x_old,fval,fval_old,gradient,torlance,iteration,iteration_max);
            if quit_flag
                done=1;
                x_best=x;
                favl_best=fval;
            end
            iteration=iteration+1;
        end
        
        % main loop for cubic interpolation
        while ~done
            
            x_base=x_old;
            x_relative=x/x_old;
            interp_matrix=[1,1,1,1;
                3,2,1,0;
                x_relative^3,x_relative^2,x_relative,1;
                3*x_relative^2,2*x_relative,1,0];
            
            if rcond(interp_matrix) < eps
                disp('error');
            end
            
            interp_value=[fval_old;gradient_old*x_base;fval;gradient*x_base];
            [x_inter_rel,coefficient_cubic]=minCubicInterpolate(interp_matrix,interp_value);
            x_inter=x_inter_rel*x_base;
            
            if DRAW_FLAG
                x_draw=1:direction*draw_interval:direction*draw_range;
                x_draw=x_draw/x_base;
                line(x_draw*x_base,coefficient_cubic(1)*x_draw.^3+coefficient_cubic(2)*x_draw.^2+...
                    coefficient_cubic(3)*x_draw+coefficient_cubic(4));
            end
            
            % limit search space, process constraints
            if x_inter > up_bou
                x_inter=up_bou;
            elseif x_inter < low_bou
                x_inter=low_bou;
            end
            
            [fval_inter,gradient_inter]=object_function(x_inter);NFE=NFE+1;
            
            % only work for one best(convex)
            % three situation discuss
            if gradient < 0
                x_old=x;
                fval_old=fval;
                gradient_old=gradient;
            else
                if gradient_inter < 0
                    x_old=x;
                    fval_old=fval;
                    gradient_old=gradient;
                end
            end
            
            x=x_inter;
            fval=fval_inter;
            gradient=gradient_inter;
            
            quit_flag=judgeQuit...
                (x,x_old,fval,fval_old,gradient,torlance,iteration,iteration_max);
            if quit_flag
                done=1;
                x_best=x;
                favl_best=fval;
            end
            
            result_x_list=[result_x_list;x];
            result_fval_list=[result_fval_list;fval];
            iteration=iteration+1;
        end
        output.result_x_list=result_x_list;
        output.result_fval_list=result_fval_list;
        
        function [lamada,coefficient_cubic]=minCubicInterpolate(interpolate_matrix,interpolate_value)
            % calculate min cubic curve
            %
            coefficient_cubic=interpolate_matrix\interpolate_value;
            
            temp_sqrt=4*coefficient_cubic(2)^2-12*coefficient_cubic(1)*coefficient_cubic(3);
            if temp_sqrt>=0
                temp_lamada=-coefficient_cubic(2)/3/coefficient_cubic(1)+...
                    sqrt(temp_sqrt)/6/coefficient_cubic(1);
                if (temp_lamada*6*coefficient_cubic(1)+2*coefficient_cubic(2))>0
                    lamada=temp_lamada;
                else
                    lamada=-coefficient_cubic(2)/3/coefficient_cubic(1)-...
                        sqrt(temp_sqrt)...
                        /6/coefficient_cubic(1);
                end
            else
                lamada=-coefficient_cubic(2)/3/coefficient_cubic(1);
            end
        end
        function quit_flag=judgeQuit...
                (x,x_old,fval,fval_old,gradient,torlance,iteration,iteration_max)
            quit_flag=0;
            if abs(fval-fval_old)/fval_old < torlance
                quit_flag=1;
            end
            if abs(gradient) < torlance
                quit_flag=1;
            end
            if abs(x-x_old) < 1e-5
                quit_flag=1;
            end
            if iteration >= iteration_max
                quit_flag=1;
            end
        end
    end
    function predict_y=interpolationEnsemleRadialBasisPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            model_number,beta_list,basis_function_list,c_list,...
            w,predict_x)
        % ensemle radial basis function interpolation predict function
        % input predict_x and radialbasis_model model
        if size(predict_x,2) > 1
            predict_x=predict_x';
        end
        
        % normalize data
        predict_x_nomlz=(predict_x-aver_X')./stdD_X';
        
        % calculate each sub model predict fval and get predict_y
        predict_y_nomlz=0;
        for model_index__=1:model_number
            basis_function__=basis_function_list{model_index__};
            c__=c_list(model_index__);
            beta__=beta_list(:,:,model_index__);
            predict_y_nomlz_sub=interpolationRadialBasisPredictor...
                (X_nomlz,beta__,basis_function__,c__,predict_x_nomlz);
            predict_y_nomlz=predict_y_nomlz+w(model_index__)*predict_y_nomlz_sub;
        end
        
        % normalize data
        predict_y=predict_y_nomlz*stdD_Y+aver_Y;
        
        function predict_y_nomlz=interpolationRadialBasisPredictor...
                (X_nomlz,beta,basis_function,c,predict_x_nomlz)
            % radial basis function interpolation predict function
            % input predict_x and radialbasis_model model
            % predict_x is row vector
            % output the predict value
            [x_number__,~]=size(X_nomlz);
            
            % predict value
            X_inner_product__=zeros(x_number__,1);
            for x_index__=1:x_number__
                X_inner_product__(x_index__)=...
                    basis_function(sum((X_nomlz(x_index__,:)'-predict_x_nomlz).^2),c);
            end
            
            % predict variance
            predict_y_nomlz=beta'*X_inner_product__;
        end
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

function drawFunction(object_function,low_bou,up_bou,...
    grid_number,Y_min,Y_max,figure_handle)
if nargin < 7
    figure_handle=figure(10);
    if nargin < 6
        Y_max=inf;
        if nargin < 5
            Y_min=-inf;
            if nargin < 4
                grid_number=100;
            end
        end
    end
end
axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end
dimension=size(low_bou,1);

switch dimension
    case 1
        d_bou=(up_bou-low_bou)/grid_number;
        X__=low_bou:d_bou:up_bou;
        fval__=zeros(grid_number+1,1);
        for x_index__=1:(grid_number+1)
            fval__(x_index__)=object_function(X__(x_index__));
        end
        line(axes_handle,X__,fval__);
        xlabel('X');
        ylabel('value');
        
    case 2
        d_bou=(up_bou-low_bou)/grid_number;
        [X__,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
        fval__=zeros(grid_number+1);
        for x_index__=1:grid_number+1
            for y_index__=1:grid_number+1
                predict_x=([x_index__;y_index__]-1).*d_bou+low_bou;
                fval__(x_index__,y_index__)=object_function(predict_x);
            end
        end
        fval__(find(fval__ > Y_max))=Y_max;
        fval__(find(fval__ < Y_min))=Y_min;
        surf(axes_handle,X__',Y',fval__,'FaceAlpha',0.5,'EdgeColor','none');
        xlabel('X');
        ylabel('Y');
        zlabel('value');
end

end
function [gradient,hessian]=differ(differ_function,x,fval,step)
% differ function to get gradient and hessian
%
variable_number=length(x);
if nargin < 4
    step=1e-6;
end
if nargin < 3
    fval=differ_function(x);
end
fval__=zeros(variable_number,2); % backward is 1, forward is 2
gradient=zeros(variable_number,1);
hessian=zeros(variable_number);

% fval and gradient
for variable_index__=1:variable_number
    x_forward__=x;
    x_backward__=x;
    x_backward__(variable_index__)=x_backward__(variable_index__)-step;
    fval__(variable_index__,1)=differ_function(x_backward__);
    
    x_forward__(variable_index__)=x_forward__(variable_index__)+step;
    fval__(variable_index__,2)=differ_function(x_forward__);
    
    gradient(variable_index__)=...
        (fval__(variable_index__,2)-fval__(variable_index__,1))/2/step;
end

% hessian
for variable_index__=1:variable_number
    hessian(variable_index__,variable_index__)=...
        (fval__(variable_index__,2)-2*fval+fval__(variable_index__,1))/step/step;
    for variable_index_next__=variable_index__+1:variable_number
        x_for_for=x;
        x_for_for(variable_index__)=x_for_for(variable_index__)+step;
        x_for_for(variable_index_next__)=x_for_for(variable_index_next__)+step;
        
        hessian(variable_index__,variable_index_next__)=(...
            differ_function(x_for_for)-...
            fval__(variable_index__,2)-fval__(variable_index_next__,2)+...
            fval...
            )/step/step;
        hessian(variable_index_next__,variable_index__)=...
            hessian(variable_index__,variable_index_next__);
    end
end
end