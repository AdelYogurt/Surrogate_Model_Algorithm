clc;
clear;
close all hidden;

% x_number=50;
% dimension=2;
% times=1;
% X1=rand(x_number,dimension)*10;
% X2=rand(x_number*times,dimension)*2+[5,5];
% class1=1*ones(x_number,1);
% class2=0*ones(x_number*times,1);
% X=[X1;X2];class=[class1;class2];

% load('matlab.mat');

func_low=@(x) lowFunction(x);
func_high=@(x) highFunction(x);

func_low_bou=@(x) 0.45+sin(2.2*pi*x)/2.5;
func_high_bou=@(x) 0.5+sin(2.5*pi*x)/3;

variable_number=2;
x_number=45;
X=zeros(x_number,2);
for x_index=1:30
    x=x_index/30-1/60;
    X(x_index,:)=[x,func_low_bou(x)+(rand()-0.5)*0.3];
end
X(31:45,:)=lhsdesign(15,variable_number);
Y=zeros(x_number,1);
for x_index=1:x_number
    Y(x_index)=func_low(X(x_index,:)');
end

x_number_H=25;
X_H=zeros(x_number_H,2);
for x_index=1:x_number_H
    x=x_index/x_number_H-1/x_number_H/2;
    X_H(x_index,:)=[x,func_high_bou(x)+(rand()-0.5)*0.3];
end
Y_H=zeros(x_number_H,1);
for x_index=1:x_number_H
    Y_H(x_index)=func_high(X_H(x_index,:)');
end

[SVM_predict_function,SVM_model]=classifySupportVectorMachine...
    (X,Y,10);
b=SVM_model.b;

classifyVisualization(SVM_model)

function [SVM_predict_function,SVM_model]=classifySupportVectorMachine...
    (X,Class,C,low_bou,up_bou,kernel_function)
% generate support vector machine model version 0
% version 0 use fmincon to get alpha
% only support binary classification, 0 and 1
% X, Y is x_number x variable_number matrix
% C is penalty factor, default is empty
% kernel_function default is gauss kernal function
% kernel_function should be @(x1,x2) ...
%
if nargin < 6
    kernel_function=[];
    if nargin < 3
        C=[];
    end
end

[x_number,variable_number]=size(X);

% normalization data
if nargin < 5
    up_bou=max(X);
    if nargin < 4
        low_bou=min(X);
    end
end
X_nomlz=(X-low_bou)./(up_bou-low_bou);

% transfer class into y
Y=Class;
for x_index=1:x_number
    if Y(x_index) == 0
        Y(x_index)=-1;
    end
end

% default kernal function
if isempty(kernel_function)
    sigma=-100*log(1/sqrt(x_number))/variable_number^2;
    kernel_function=@(U,V) kernelGaussian(U,V,sigma);
end

% initialization kernal function process X_cov
K=kernel_function(X_nomlz,X_nomlz);

% min SVM object function to get alpha
object_function_SVM=@(alpha) -objectFunctionSVM(alpha,K,Y);
alpha_initial=ones(x_number,1)*0.5;
low_bou_fmincon=0*ones(x_number,1);
if isempty(C) || C==0
    up_bou_fmincon=[];
else
    up_bou_fmincon=C*ones(x_number,1);
end
Aeq=Y';
fmincon_options = optimoptions('fmincon','Display','none','Algorithm','sqp');
alpha=fmincon(object_function_SVM,alpha_initial,...
    [],[],Aeq,0,low_bou_fmincon,up_bou_fmincon,[],fmincon_options);

% obtain other paramter
alpha_Y=alpha.*Y;

w=sum(alpha_Y.*X_nomlz);
index_list=find(alpha > 1e-6); % support vector
alpha_Y_cov=K*alpha_Y;
b=sum(Y(index_list)-alpha_Y_cov(index_list))/length(index_list);

% generate predict function
SVM_predict_function=@(x) classifySupportVectorMachinePredictor...
    (x,X_nomlz,alpha_Y,b,low_bou,up_bou,kernel_function);

% output model
SVM_model.X=X;
SVM_model.Class=Class;
SVM_model.Y=Y;
SVM_model.X_nomlz=X_nomlz;
SVM_model.low_bou=low_bou;
SVM_model.up_bou=up_bou;
SVM_model.alpha=alpha;
SVM_model.w=w;
SVM_model.b=b;
SVM_model.kernel_function=kernel_function;
SVM_model.predict_function=SVM_predict_function;

    function fval=objectFunctionSVM(alpha,X_inner_product,Y)
        % support vector machine maximum object function
        %
        alpha=alpha(:);
        alpha_Y__=alpha.*Y;
        fval=sum(alpha)-alpha_Y__'*X_inner_product*alpha_Y__/2;
    end
    function [Class_pred,Probability]=classifySupportVectorMachinePredictor...
            (X_pred,X_nomlz,alpha_Y,b,low_bou,up_bou,kernel_function)
        % predict_fval is 1 or -1, predict_class is 1 or 0
        %
        % x input is colume vector
        %
        X_pred_nomlz=(X_pred-low_bou)./(up_bou-low_bou);
        K_pred=kernel_function(X_pred_nomlz,X_nomlz);
        Probability=K_pred*alpha_Y+b;
        Probability=1./(1+exp(-Probability));
        Class_pred=Probability > 0.5;
    end
    function K=kernelGaussian(U,V,sigma)
        % gaussian kernal function
        %
        K=zeros(size(U,1),size(V,1));
        vari_num=size(U,2);
        for vari_index=1:vari_num
            K=K+(U(:,vari_index)-V(:,vari_index)').^2;
        end
        K=exp(-K*sigma);
    end
end

function fval=lowFunction(x)
if 0.45+sin(2.2*pi*x(1))/2.5-x(2) > 0
    fval=1;
else
    fval=-1;
end
end
function fval=highFunction(x)
if 0.5+sin(2.5*pi*x(1))/3-x(2) > 0
    fval=1;
else
    fval=-1;
end
end