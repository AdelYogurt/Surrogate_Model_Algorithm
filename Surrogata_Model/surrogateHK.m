clc;
% clear;
close all hidden;

% load('MFK.mat')
% 
% K_fold=5;
% 
% HF_select=randi(100,1,size(matPntHF,1)/K_fold*(K_fold-1));
% HF_check=1:100;
% HF_check(HF_select)=[];
% LF_select=randi(400,1,size(matPntLF1,1)/K_fold*(K_fold-1));
% 
% XHF=matPntHF(HF_select,:);
% YHF=matObjHF(HF_select,1);
% 
% XLF=matPntLF1(LF_select,:);
% YLF=matObjLF1(LF_select,1);
% 
% [predict_function,MK_model]=interpMultiKrigingPreModel...
%     (XHF, YHF, XLF, YLF);
% 
% check_error_list=zeros(size(matPntHF,1)/K_fold,1);
% for check_index=1:(size(matPntHF,1)/K_fold)
%     x_index=HF_check(check_index);
%     check_error_list(check_index)=(predict_function(matPntHF(x_index,:))-matObjHF(x_index,1));
% end
% disp(['max error: ',num2str(max(check_error_list))]);

benchmark=BenchmarkFunction();
xLF_number=11;
XHF=[0;0.4;0.6;1];
YHF=benchmark.singleForresterObject(XHF);
XLF=linspace(0,1,xLF_number)';
YLF=benchmark.singleForresterObjectLow(XLF);
[predict_function_HK,HK_model]=interpHieraKrigingPreModel...
    (XHF,YHF,XLF,YLF);
[predict_function,HK_model]=interpKrigingPreModel(XHF,YHF);
x=(0:0.01:1)';
[y_HK]=predict_function_HK(x);
[y]=predict_function(x);
line(x,benchmark.singleForresterObject(x),'Color','k','LineStyle','-')
line(x,benchmark.singleForresterObjectLow(x),'Color','k','LineStyle','--');
line(x,y_HK,'Color','b','LineStyle','-')
line(x,y,'Color','b','LineStyle','--');
legend({'singleForresterObject','singleForresterObjectLow','HK','K'})

function [predict_function,HK_model]=interpHieraKrigingPreModel...
    (XHF,YHF,XLF,YLF)
% construct Hierarchical Kriging version 0
% XHF, YHF are x_HF_number x variable_number matrix
% XLF, YLF are x_LF_number x variable_number matrix
% aver_X,stdD_X is 1 x x_HF_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
%
% input:
% XHF, YHF, XLF, YLF
%
% output:
% kriging model (predict_function,
% XHF, YHF, XLF, YLF, base_function_list)
%
% reference: [1] HAN Z-H, GöRTZ S. Hierarchical Kriging Model for
% Variable-Fidelity Surrogate Modeling [J]. AIAA Journal, 2012, 50(9):
% 1885-96.
%
% Copyright 2023.2 Adel
%
X=[XHF;XLF];
Y=[YHF;YLF];
[x_number,variable_number]=size(X);
x_HF_number=size(XHF,1);
x_LF_number=size(XLF,1);
if nargin < 5
    theta=ones(1,variable_number);
end

% normalize data
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
index__=find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__)=1; end
index__=find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__)=1; end
XHF_nomlz=(XHF-repmat(aver_X,x_HF_number,1))./repmat(stdD_X,x_HF_number,1);
YHF_nomlz=(YHF-repmat(aver_Y,x_HF_number,1))./repmat(stdD_Y,x_HF_number,1);
XLF_nomlz=(XLF-repmat(aver_X,x_LF_number,1))./repmat(stdD_X,x_LF_number,1);
YLF_nomlz=(YLF-repmat(aver_Y,x_LF_number,1))./repmat(stdD_Y,x_LF_number,1);

% first step
% construct low fidelity model

% initial X_dis_sq
XHF_dis_sq=zeros(x_HF_number,x_HF_number,variable_number);
for variable_index=1:variable_number
    XHF_dis_sq(:,:,variable_index)=...
        (XHF_nomlz(:,variable_index)-XHF_nomlz(:,variable_index)').^2;
end
XLF_dis_sq=zeros(x_LF_number,x_LF_number,variable_number);
for variable_index=1:variable_number
    XLF_dis_sq(:,:,variable_index)=...
        (XLF_nomlz(:,variable_index)-XLF_nomlz(:,variable_index)').^2;
end

% regression function define
% reg_function=@(X) regZero(X);
reg_function=@(X) regLinear(X);

% calculate reg
fval_reg_LF=reg_function(XLF_nomlz);

% optimal to get hyperparameter
fmincon_option=optimoptions('fmincon','Display','iter-detailed',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',true);
low_bou_kriging=1e-1*ones(variable_number,1);
up_bou_kriging=20*ones(variable_number,1);
object_function=@(theta) objectFunctionKriging...
    (XLF_dis_sq,YLF_nomlz,x_LF_number,variable_number,theta,fval_reg_LF);

theta=fmincon...
    (object_function,theta,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% get parameter
[cov_LF,inv_cov_LF,beta_LF,sigma_sq_LF]=interpKriging...
    (XLF_dis_sq,YLF_nomlz,x_LF_number,variable_number,theta,fval_reg_LF);
gama_LF=inv_cov_LF*(YLF_nomlz-fval_reg_LF*beta_LF);
FTRF_LF=fval_reg_LF'*inv_cov_LF*fval_reg_LF;

% initialization predict function
predict_function_LF=@(X_predict) interpKrigingPredictor...
    (XLF_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_LF_number,variable_number,theta,beta_LF,gama_LF,sigma_sq_LF,...
    inv_cov_LF,fval_reg_LF,FTRF_LF,X_predict,reg_function);

% second step
% construct hierarchical model

% evaluate low fidelty predict value in high fidelity point as base fval
reg_function=@(X) predict_function_LF(X);

% calculate reg
fval_reg_HF=reg_function(XHF);

% optimal to get hyperparameter
object_function=@(theta) objectFunctionKriging...
    (XHF_dis_sq,YHF_nomlz,x_HF_number,variable_number,theta,fval_reg_HF);

% [fval,gradient]=object_function(theta)
% [~,gradient_differ]=differ(object_function,theta)

theta=fmincon...
    (object_function,theta,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% calculate covariance and other parameter
[cov_HF,inv_cov_HF,beta_HF,sigma_sq_HF]=interpKriging...
    (XHF_dis_sq,YHF_nomlz,x_HF_number,variable_number,theta,fval_reg_HF);
gama_HF=inv_cov_HF*(YHF_nomlz-fval_reg_HF*beta_HF);
FTRF_HF=fval_reg_HF'*inv_cov_HF*fval_reg_HF;

% initialization predict function
predict_function=@(X_predict) interpKrigingPredictor...
    (XHF_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_HF_number,variable_number,theta,beta_HF,gama_HF,sigma_sq_HF,...
    inv_cov_HF,fval_reg_HF,FTRF_HF,X_predict,reg_function);

HK_model.X=XHF;
HK_model.Y=YHF;
HK_model.XHF=XHF;
HK_model.YHF=YHF;
HK_model.XLF=XLF;
HK_model.YLF=YLF;
HK_model.cov_LF=cov_LF;
HK_model.inv_cov_LF=inv_cov_LF;
HK_model.cov_HF=cov_HF;
HK_model.inv_cov_HF=inv_cov_HF;

HK_model.theta=theta;
HK_model.aver_X=aver_X;
HK_model.stdD_X=stdD_X;
HK_model.aver_Y=aver_Y;
HK_model.stdD_Y=stdD_Y;

HK_model.predict_function=predict_function;

% abbreviation:
% num: number, pred: predict, vari: variable
    function [fval,gradient]=objectFunctionKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg)
        % function to minimize sigma_sq
        %
        [cov,inv_cov,~,sigma2,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg);
        
        % calculation negative log likelihood
        L=chol(cov)';
        fval=x_num/2*log(sigma2)+sum(log(diag(L)));

        % calculate gradient
        if nargout > 1
            % gradient
            gradient=zeros(vari_num,1);
            for vari_index=1:vari_num
                dcov_dtheta=-(X_dis_sq(:,:,vari_index).*cov);

                dinv_cov_dtheta=...
                    -inv_cov*dcov_dtheta*inv_cov;

                dinv_FTRF_dtheta=-inv_FTRF*...
                    (F_reg'*dinv_cov_dtheta*F_reg)*...
                    inv_FTRF;
                
                dmiu_dtheta=dinv_FTRF_dtheta*(F_reg'*inv_cov*Y)+...
                    inv_FTRF*(F_reg'*dinv_cov_dtheta*Y);
                
                dY_Fmiu_dtheta=-F_reg*dmiu_dtheta;

                dsigma2_dtheta=(dY_Fmiu_dtheta'*inv_cov*Y_Fmiu+...
                    Y_Fmiu'*dinv_cov_dtheta*Y_Fmiu+...
                    Y_Fmiu'*inv_cov*dY_Fmiu_dtheta)/x_num;
                
                dlnsigma2_dtheta=1/sigma2*dsigma2_dtheta;

                dlndetR=trace(inv_cov*dcov_dtheta);

                gradient(vari_index)=x_num/2*dlnsigma2_dtheta+0.5*dlndetR;
            end
        end
    end
    function [cov,inv_cov,beta,sigma_sq,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg)
        % kriging interpolation kernel function
        % Y(x)=beta+Z(x)
        %
        cov=zeros(x_num,x_num);
        for vari_index=1:vari_num
            cov=cov+X_dis_sq(:,:,vari_index)*theta(vari_index);
        end
        cov=exp(-cov)+eye(x_num)*1e-6;

        % coefficient calculation
        inv_cov=cov\eye(x_num);
        inv_FTRF=(F_reg'*inv_cov*F_reg)\eye(size(F_reg,2));

        % basical bias
        beta=inv_FTRF*(F_reg'*inv_cov*Y);
        Y_Fmiu=Y-F_reg*beta;
        sigma_sq=(Y_Fmiu'*inv_cov*Y_Fmiu)/x_num;
        
    end
    function [Y_pred,Var_pred]=interpKrigingPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,theta,beta,gama,sigma_sq,...
            inv_cov,fval_reg,FTRF,X_pred,reg_function)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num,~]=size(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;
        
        % predict covariance
        predict_cov=zeros(x_num,x_pred_num);
        for vari_index=1:vari_num
            predict_cov=predict_cov+...
                (X_nomlz(:,vari_index)-X_pred_nomlz(:,vari_index)').^2*theta(vari_index);
        end
        predict_cov=exp(-predict_cov);

        % predict base fval
        predict_fval_reg=reg_function(X_pred);
        Y_pred=predict_fval_reg*beta+predict_cov'*gama;
        
        % predict variance
        u__=fval_reg'*inv_cov*predict_cov-predict_fval_reg';
        Var_pred=sigma_sq*...
            (1+u__'/FTRF*u__+...
            -predict_cov'*inv_cov*predict_cov);
        
        % normalize data
        Y_pred=Y_pred*stdD_Y+aver_Y;
        Var_pred=diag(Var_pred)*stdD_Y*stdD_Y;
    end
    function F_reg=regZero(X)
        % zero order base funcion
        %
        F_reg=ones(size(X,1),1); % zero
    end
    function F_reg=regLinear(X)
        % first order base funcion
        %
        F_reg=[ones(size(X,1),1),X]; % linear
    end
end

function [predict_function,kriging_model]=interpKrigingPreModel...
    (X,Y,theta)
% version 6, nomalization method is grassian
% add multi x_predict input support
% improve constrcut speed
% prepare model, optimal theta and calculation parameter
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
%
% input initial data X, Y, which are real data
%
% output is a kriging model, include predict_function...
% X, Y, base_function_list
%
% Copyright 2023.2 Adel
%
[x_number,variable_number]=size(X);
if nargin < 3
    theta=ones(1,variable_number);
end

% normalize data
aver_X=mean(X);
stdD_X=std(X);
aver_Y=mean(Y);
stdD_Y=std(Y);
index__=find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__)=1; end
index__=find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__)=1; end
X_nomlz=(X-aver_X)./stdD_X;
Y_nomlz=(Y-aver_Y)./stdD_Y;

% initial X_dis_sq
X_dis_sq=zeros(x_number,x_number,variable_number);
for variable_index=1:variable_number
    X_dis_sq(:,:,variable_index)=...
        (X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end

% regression function define
% reg_function=@(X) regZero(X);
reg_function=@(X) regLinear(X);

% calculate reg
fval_reg=reg_function(X_nomlz);

% optimal to get hyperparameter
fmincon_option=optimoptions('fmincon','Display','iter-detailed',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',true);
low_bou_kriging=1e-1*ones(variable_number,1);
up_bou_kriging=20*ones(variable_number,1);
object_function_kriging=@(theta) objectFunctionKriging...
    (X_dis_sq,Y_nomlz,x_number,variable_number,theta,fval_reg);

% [fval,gradient]=object_function_kriging(theta)
% [~,gradient_differ]=differ(object_function_kriging,theta)

theta=fmincon...
    (object_function_kriging,theta,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% get parameter
[covariance,inv_covariance,beta,sigma_sq]=interpKriging...
    (X_dis_sq,Y_nomlz,x_number,variable_number,theta,fval_reg);
gama=inv_covariance*(Y_nomlz-fval_reg*beta);
FTRF=fval_reg'*inv_covariance*fval_reg;

% initialization predict function
predict_function=@(predict_x) interpKrigingPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,theta,beta,gama,sigma_sq,...
    inv_covariance,fval_reg,FTRF,predict_x,reg_function);

kriging_model.X=X;
kriging_model.Y=Y;
kriging_model.X_normalize=X_nomlz;
kriging_model.Y_normalize=Y_nomlz;
kriging_model.fval_regression=fval_reg;
kriging_model.covariance=covariance;
kriging_model.inv_covariance=inv_covariance;

kriging_model.theta=theta;
kriging_model.beta=beta;
kriging_model.gama=gama;
kriging_model.sigma_sq=sigma_sq;
kriging_model.aver_X=aver_X;
kriging_model.stdD_X=stdD_X;
kriging_model.aver_Y=aver_Y;
kriging_model.stdD_Y=stdD_Y;

kriging_model.predict_function=predict_function;

% abbreviation:
% num: number, pred: predict, vari: variable
    function [fval,gradient]=objectFunctionKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg)
        % function to minimize sigma_sq
        %
        [cov,inv_cov,~,sigma2,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg);
        
%         % calculation negative log likelihood
%         L=chol(cov)';
%         fval_reg_L=L\fval_reg__;
%         beta__=(fval_reg_L'*fval_reg_L)\(fval_reg_L'*(L\Y));
%         Y_fbeta=Y-fval_reg__*beta__;
%         Y_fbeta_L=L\Y_fbeta;
% %         fval=(Y_fbeta_L'*Y_fbeta_L)/x_num;
%         fval=x_num/2*log((Y_fbeta_L'*Y_fbeta_L)/x_num)+0.5*sum(log(diag(L)));

        % calculation negative log likelihood
        L=chol(cov)';
        fval=x_num/2*log(sigma2)+sum(log(diag(L)));

        % calculate gradient
        if nargout > 1
            % gradient
            gradient=zeros(vari_num,1);
            for vari_index=1:vari_num
                dcov_dtheta=-(X_dis_sq(:,:,vari_index).*cov);

                dinv_cov_dtheta=...
                    -inv_cov*dcov_dtheta*inv_cov;

                dinv_FTRF_dtheta=-inv_FTRF*...
                    (F_reg'*dinv_cov_dtheta*F_reg)*...
                    inv_FTRF;
                
                dmiu_dtheta=dinv_FTRF_dtheta*(F_reg'*inv_cov*Y)+...
                    inv_FTRF*(F_reg'*dinv_cov_dtheta*Y);
                
                dY_Fmiu_dtheta=-F_reg*dmiu_dtheta;

                dsigma2_dtheta=(dY_Fmiu_dtheta'*inv_cov*Y_Fmiu+...
                    Y_Fmiu'*dinv_cov_dtheta*Y_Fmiu+...
                    Y_Fmiu'*inv_cov*dY_Fmiu_dtheta)/x_num;
                
                dlnsigma2_dtheta=1/sigma2*dsigma2_dtheta;

                dlndetR=trace(inv_cov*dcov_dtheta);

                gradient(vari_index)=x_num/2*dlnsigma2_dtheta+0.5*dlndetR;
            end
        end
    end
    function [cov,inv_cov,beta,sigma_sq,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,theta,F_reg)
        % kriging interpolation kernel function
        % Y(x)=beta+Z(x)
        %
        cov=zeros(x_num,x_num);
        for vari_index=1:vari_num
            cov=cov+X_dis_sq(:,:,vari_index)*theta(vari_index);
        end
        cov=exp(-cov)+eye(x_num)*1e-6;

        % coefficient calculation
        inv_cov=cov\eye(x_num);
        inv_FTRF=(F_reg'*inv_cov*F_reg)\eye(size(F_reg,2));

        % basical bias
        beta=inv_FTRF*(F_reg'*inv_cov*Y);
        Y_Fmiu=Y-F_reg*beta;
        sigma_sq=(Y_Fmiu'*inv_cov*Y_Fmiu)/x_num;
        
    end
    function [Y_pred,Var_pred]=interpKrigingPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,theta,beta,gama,sigma_sq,...
            inv_cov,fval_reg,FTRF,X_pred,reg_function)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num,~]=size(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;
        
        % predict covariance
        predict_cov=zeros(x_num,x_pred_num);
        for vari_index=1:vari_num
            predict_cov=predict_cov+...
                (X_nomlz(:,vari_index)-X_pred_nomlz(:,vari_index)').^2*theta(vari_index);
        end
        predict_cov=exp(-predict_cov);

        % predict base fval
        predict_fval_reg=reg_function(X_pred_nomlz);
        Y_pred=predict_fval_reg*beta+predict_cov'*gama;
        
        % predict variance
        u__=fval_reg'*inv_cov*predict_cov-predict_fval_reg';
        Var_pred=sigma_sq*...
            (1+u__'/FTRF*u__+...
            -predict_cov'*inv_cov*predict_cov);
        
        % normalize data
        Y_pred=Y_pred*stdD_Y+aver_Y;
        Var_pred=diag(Var_pred)*stdD_Y*stdD_Y;
    end
    function F_reg=regZero(X)
        % zero order base funcion
        %
        F_reg=ones(size(X,1),1); % zero
    end
    function F_reg=regLinear(X)
        % first order base funcion
        %
        F_reg=[ones(size(X,1),1),X]; % linear
    end
end
