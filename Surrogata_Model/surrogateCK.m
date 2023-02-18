clc;
clear;
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
[predict_function_MK,MK_model]=interpCoKrigingPreModel...
    (XHF, YHF, XLF, YLF);
[predict_function,K_model]=interpKrigingPreModel(XHF, YHF);
x=(0:0.01:1)';
[y_MK]=predict_function_MK(x);
[y]=predict_function(x);
line(x,benchmark.singleForresterObject(x),'Color','k','LineStyle','-')
line(x,benchmark.singleForresterObjectLow(x),'Color','k','LineStyle','--');
line(x,y_MK,'Color','b','LineStyle','-')
line(x,y,'Color','b','LineStyle','--');
legend({'singleForresterObject','singleForresterObjectLow','MK','K'})

function [predict_function,MK_model]=interpCoKrigingPreModel...
    (XHF,YHF,XLF,YLF,hyp)
% construct Co-Kriging version 0
% XHF, YHF are x_HF_number x variable_number matrix
% XLF, YLF are x_LF_number x variable_number matrix
% aver_X,stdD_X is 1 x x_HF_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
%
% input:
% XHF, YHF, XLF, YLF, hyp(theta_HF, rho, theta_LF)
%
% output:
% kriging model (predict_function,
% XHF, YHF, XLF, YLF)
%
% reference: [1] FORRESTER A I J, SóBESTER A, KEANE A J. Multi-fidelity
% optimization via surrogate modelling [J]. Proceedings of the Royal
% Society A: Mathematical, Physical and Engineering Sciences, 2007,
% 463(3251 - 69.
%
% Copyright 2023.2 Adel
%
X=[XHF;XLF];
Y=[YHF;YLF];
[x_number,variable_number]=size(X);
x_HF_number=size(XHF,1);
x_LF_number=size(XLF,1);
if nargin < 5
    hyp=ones(1,2*variable_number+1);
end
hyp_LF=hyp((variable_number+2):end);
hyp_d=hyp(1:variable_number+1);

% normalize data
aver_X=mean(X);
stdD_X=std(X);
aver_Y=mean(Y);
stdD_Y=std(Y);
% aver_X=0;
% stdD_X=1;
% aver_Y=0;
% stdD_Y=1;
index__=find(stdD_X == 0);
if  ~isempty(index__),  stdD_X(index__)=1; end
index__=find(stdD_Y == 0);
if  ~isempty(index__),  stdD_Y(index__)=1; end
XHF_nomlz=(XHF-aver_X)./stdD_X;
YHF_nomlz=(YHF-aver_Y)./stdD_Y;
XLF_nomlz=(XLF-aver_X)./stdD_X;
YLF_nomlz=(YLF-aver_Y)./stdD_Y;
X_nomlz=[XHF_nomlz;XLF_nomlz];
Y_nomlz=[YHF_nomlz;YLF_nomlz];

% first step
% construct low fidelity model

% initial X_dis_sq
X_dis_sq=zeros(x_number,x_number,variable_number);
for variable_index=1:variable_number
    X_dis_sq(:,:,variable_index)=...
        (X_nomlz(:,variable_index)-X_nomlz(:,variable_index)').^2;
end
XHF_dis_sq=X_dis_sq(1:x_HF_number,1:x_HF_number,:);
XLF_dis_sq=X_dis_sq(x_HF_number+1:end,x_HF_number+1:end,:);

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
low_bou_kriging=1e-1*ones(1,variable_number);
up_bou_kriging=20*ones(1,variable_number);
object_function=@(hyp) objectFunctionKriging...
    (XLF_dis_sq,YLF_nomlz,x_LF_number,variable_number,fval_reg_LF,hyp);

hyp_LF=fmincon...
    (object_function,hyp_LF,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% get low fidelity parameter
[cov_LF,inv_cov_LF,beta_LF,sigma_sq_LF]=interpKriging...
    (XLF_dis_sq,YLF_nomlz,x_LF_number,variable_number,hyp_LF,fval_reg_LF);
gama_LF=inv_cov_LF*(YLF_nomlz-fval_reg_LF*beta_LF);
FTRF_LF=fval_reg_LF'*inv_cov_LF*fval_reg_LF;

% initialization predict function
predict_function_LF=@(X_predict) interpKrigingPredictor...
    (X_predict,XLF_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_LF_number,variable_number,hyp,beta_LF,gama_LF,sigma_sq_LF,...
    inv_cov_LF,fval_reg_LF,FTRF_LF,reg_function);

% second step
% construct bias model
% notice rho is hyperparamter

% evaluate error in high fidelity point
YHF_pred=zeros(x_HF_number,1); % 
for x_index=1:x_HF_number
    YHF_pred(x_index,1)=...
        predict_function_LF(XHF(x_index,:));
end
YHF_pred_nomlz=(YHF_pred-aver_Y)./stdD_Y;

% calculate reg
fval_reg_HF=reg_function(XHF_nomlz);

% optimal to get hyperparameter
fmincon_option=optimoptions('fmincon','Display','iter-detailed',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',true);
low_bou_kriging=1e-1*ones(1,variable_number+1);
up_bou_kriging=20*ones(1,variable_number+1);
object_function=@(hyp) objectFunctionCoKriging...
    (XHF_dis_sq,YHF_nomlz,YHF_pred_nomlz,x_HF_number,variable_number,fval_reg_HF,hyp);

% [fval,gradient]=object_function(hyp_HF)
% [~,gradient_differ]=differ(object_function,hyp_HF)

hyp_d=fmincon...
    (object_function,hyp_d,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% get high fidelity parameter
D_nomlz=YHF_nomlz-hyp_d(end)*YHF_pred_nomlz;
[cov_d,inv_cov_d,beta_d,sigma_sq_d]=interpKriging...
    (XHF_dis_sq,D_nomlz,x_HF_number,variable_number,hyp_d,fval_reg_HF);
% gama_HF=inv_cov_HF*(D_nomlz-fval_reg_HF*beta_HF);
% FTRF_HF=fval_reg_HF'*inv_cov_HF*fval_reg_HF;

% get total model parameter
covariance=calCoCov(X_dis_sq,x_number,variable_number,[hyp_d,hyp_LF],sigma_sq_d,sigma_sq_LF);
inv_covariance=covariance\eye(x_number);
fval_reg=reg_function(X_nomlz);
beta=(fval_reg'*inv_covariance*fval_reg)\(fval_reg'*inv_covariance*Y_nomlz);
gama=inv_covariance*(Y_nomlz-fval_reg*beta);

% initialization predict function
predict_function=@(X_pred) interpCoKrigingPredictor...
    (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,x_HF_number,hyp_d(1:variable_number),hyp_d(end),hyp_LF,beta,gama,sigma_sq_d,sigma_sq_LF,...
    inv_covariance,reg_function);

MK_model.X=XHF;
MK_model.Y=YHF;
MK_model.XHF=XHF;
MK_model.YHF=YHF;
MK_model.XLF=XLF;
MK_model.YLF=YLF;
MK_model.predict_function=predict_function;


% abbreviation:
% num: number, pred: predict, vari: variable
    function [fval,gradient]=objectFunctionKriging...
            (X_dis_sq,Y,x_num,vari_num,F_reg,hyp)
        % function to minimize sigma_sq
        %
        [cov,inv_cov,~,sigma2,inv_FTRF,Y_Fmiu]=interpKriging...
            (X_dis_sq,Y,x_num,vari_num,hyp,F_reg);

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
                
                dbeta_dtheta=dinv_FTRF_dtheta*(F_reg'*inv_cov*Y)+...
                    inv_FTRF*(F_reg'*dinv_cov_dtheta*Y);
                
                dY_Fmiu_dtheta=-F_reg*dbeta_dtheta;

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
    function [fval,gradient]=objectFunctionCoKriging...
            (XHF_dis_sq,YHF,YHF_pred,x_num,vari_num,F_reg,hyp)
        % function to minimize negative likelihood
        %
        theta=hyp(1:vari_num);
        rho=hyp(end);
        d=YHF-rho*YHF_pred; % notice Y equal to d

        [cov,inv_cov,~,sigma2,inv_FTRF,Y_Fmiu]=interpKriging...
            (XHF_dis_sq,d,x_num,vari_num,theta,F_reg);

        % calculation negative log likelihood
        L=chol(cov)';
        fval=x_num/2*log(sigma2)+sum(log(diag(L)));

        % calculate gradient
        if nargout > 1
            % gradient
            gradient=zeros(vari_num+1,1);

            % theta H 
            for vari_index=1:vari_num
                dcov_dtheta=-(XHF_dis_sq(:,:,vari_index).*cov);

                dinv_cov_dtheta=...
                    -inv_cov*dcov_dtheta*inv_cov;

                dinv_FTRF_dtheta=-inv_FTRF*...
                    (F_reg'*dinv_cov_dtheta*F_reg)*...
                    inv_FTRF;
                
                dbeta_dtheta=dinv_FTRF_dtheta*(F_reg'*inv_cov*YHF)+...
                    inv_FTRF*(F_reg'*dinv_cov_dtheta*YHF);
                
                dY_Fmiu_dtheta=-F_reg*dbeta_dtheta;

                dsigma2_dtheta=(dY_Fmiu_dtheta'*inv_cov*Y_Fmiu+...
                    Y_Fmiu'*dinv_cov_dtheta*Y_Fmiu+...
                    Y_Fmiu'*inv_cov*dY_Fmiu_dtheta)/x_num;
                
                dlnsigma2_dtheta=1/sigma2*dsigma2_dtheta;

                dlndetR=trace(inv_cov*dcov_dtheta);

                gradient(vari_index)=x_num/2*dlnsigma2_dtheta+0.5*dlndetR;
            end

            % rho
            dY_drho=-YHF_pred;

            dbeta_drho=inv_FTRF*(F_reg'*inv_cov*dY_drho);

            dY_Fmiu_drho=(dY_drho-F_reg*dbeta_drho);

            dsigma2_drho=(dY_Fmiu_drho'*inv_cov*Y_Fmiu+...
                Y_Fmiu'*inv_cov*dY_Fmiu_drho)/x_num;

            dlnsigma2_drho=1/sigma2*dsigma2_drho;

            gradient(end)=x_num/2*dlnsigma2_drho;
        end
    end

    function [cov,dcov_dtheta]=calCoCov(X_dis_sq,x_num,vari_num,hyp,sigma2_HF,sigma2_LF)
        % calculate covariance of x with multi fidelity
        % hyp: theta_H, theta_L, rho
        %
        theta_HF=hyp(1:vari_num);
        rho=hyp(vari_num+1);
        theta_LF=hyp((vari_num+2):end);
        
        % exp of x__x with theta H
        cov_H=zeros(x_HF_number);
        for vari_index=1:vari_num
            cov_H=cov_H+...
                X_dis_sq(1:x_HF_number,1:x_HF_number,vari_index)*theta_HF(vari_index);
        end
        cov_H=sigma2_HF*exp(-cov_H);

        % exp of x__x with theta L
        exp_disL=zeros(x_num);
        for vari_index=1:vari_num
            exp_disL=exp_disL+...
                X_dis_sq(:,:,vari_index)*theta_LF(vari_index);
        end
        exp_disL=sigma2_LF*exp(-exp_disL);
        % times rho: HH to rho2, HL to rho, LL to 1
        rho_exp_disL=exp_disL;
        rho_exp_disL(1:x_HF_number,1:x_HF_number)=...
            (rho*rho)*exp_disL(1:x_HF_number,1:x_HF_number);
        rho_exp_disL((x_HF_number+1):end,1:end)=...
            rho*exp_disL((x_HF_number+1):end,1:end);
        rho_exp_disL(1:end,(x_HF_number+1):end)=...
            rho_exp_disL((x_HF_number+1):end,1:end)';

        cov_L=rho_exp_disL;
        cov=cov_L;
        cov(1:x_HF_number,1:x_HF_number)=cov(1:x_HF_number,1:x_HF_number)+cov_H;

    end

    function [Y_pred,Var_pred]=interpKrigingPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,theta,beta,gama,sigma_sq,...
            inv_cov,fval_reg,FTRF,reg_function)
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

    function [Y_pred,Var_pred]=interpCoKrigingPredictor...
            (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,x_HF_num,theta_d,rho,theta_LF,beta,gama,sigma_sq_d,sigma_sq_LF,...
            inv_cov,reg_function)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num,~]=size(X_pred);

        % normalize data
        X_pred_nomlz=(X_pred-aver_X)./stdD_X;
        
        % distance sq of X_predict
        X_pred_dis_sq=zeros(x_num,x_pred_num,vari_num);
        for vari_index=1:vari_num
            X_pred_dis_sq(:,:,vari_index)=...
                (X_nomlz(:,vari_index)-X_pred_nomlz(:,vari_index)').^2;
        end

        % delta
        exp_dis_d=zeros(x_HF_num,x_pred_num);
        for vari_index=1:vari_num
            exp_dis_d=exp_dis_d+...
                X_pred_dis_sq(1:x_HF_num,:,vari_index)*theta_d(vari_index);
        end
        exp_dis_d=exp(-exp_dis_d);

        % LF
        exp_dis_LF=zeros(x_num,x_pred_num);
        for vari_index=1:vari_num
            exp_dis_LF=exp_dis_LF+...
                X_pred_dis_sq(:,:,vari_index)*theta_LF(vari_index);
        end
        exp_dis_LF=exp(-exp_dis_LF);

        % covariance of X_predict
        predict_cov=exp_dis_LF*rho*sigma_sq_LF;
        predict_cov(1:x_HF_num,:)=predict_cov(1:x_HF_num,:)*rho+...
            sigma_sq_d*exp_dis_d;

        % predict base fval
        predict_fval_reg=reg_function(X_pred_nomlz);
        Y_pred=predict_fval_reg*beta+predict_cov'*gama;
        
        % predict variance
        Var_pred=sigma_sq_LF*rho*rho+sigma_sq_d+...
            -predict_cov'*inv_cov*predict_cov;
        
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
X_nomlz=(X-repmat(aver_X,x_number,1))./repmat(stdD_X,x_number,1);
Y_nomlz=(Y-repmat(aver_Y,x_number,1))./repmat(stdD_Y,x_number,1);

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
