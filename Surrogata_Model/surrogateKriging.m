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
% 
% X=matPntHF(HF_select,:);
% Y=matObjHF(HF_select,1);
% 
% [predict_function,kriging_model]=interpKrigingPreModel(X,Y);
% check_error_list=zeros(size(matPntHF,1)/K_fold,1);
% for check_index=1:(size(matPntHF,1)/K_fold)
%     x_index=HF_check(check_index);
%     check_error_list(check_index)=(predict_function(matPntHF(x_index,:))-matObjHF(x_index,1));
% end
% disp(['max error: ',num2str(max(check_error_list))]);

% load('PK.mat');
% [predict_function,kriging_model]=interpKrigingPreModel(X,Y);
% figure_handle=figure(1);
% interpVisualize(kriging_model,low_bou,up_bou,figure_handle)

benchmark=BenchmarkFunction();
low_bou=[-2,-2];
up_bou=[2,2];
X=lhsdesign(10,2).*(up_bou-low_bou)+low_bou;
Y=benchmark.singleGPObject(X);
[predict_function,kriging_model]=interpKrigingPreModel(X,Y);
figure_handle=figure(1);
interpVisualize(kriging_model,low_bou,up_bou,figure_handle)

function [predict_function,kriging_model]=interpKrigingPreModel...
    (X,Y,theta)
% version 5, nomalization method is grassian
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
aver_X = mean(X);
stdD_X = std(X);
aver_Y = mean(Y);
stdD_Y = std(Y);
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

% optimal to get hyperparameter
fmincon_option=optimoptions('fmincon','Display','iter-detailed',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',true);
low_bou_kriging=1e-1*ones(variable_number,1);
up_bou_kriging=20*ones(variable_number,1);
object_function_kriging=@(theta) objectFunctionKriging...
    (X_dis_sq,X_nomlz,Y_nomlz,x_number,variable_number,theta);

% [fval,gradient]=object_function_kriging(theta)
% [~,gradient_differ]=differ(object_function_kriging,theta)

theta=fmincon...
    (object_function_kriging,theta,[],[],[],[],low_bou_kriging,up_bou_kriging,[],fmincon_option);

% get parameter
[covariance,inv_covariance,fval_reg,beta,sigma_sq]=interpKriging...
    (X_dis_sq,X_nomlz,Y_nomlz,x_number,variable_number,theta);
gama=inv_covariance*(Y_nomlz-fval_reg*beta);
FTRF=fval_reg'*inv_covariance*fval_reg;

% initialization predict function
predict_function=@(predict_x) interpKrigingPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    x_number,variable_number,theta,beta,gama,sigma_sq,...
    inv_covariance,fval_reg,FTRF,predict_x);

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
            (X_dis_sq,X,Y,x_num,vari_num,theta)
        % function to minimize sigma_sq
        %
        cov=zeros(x_num,x_num);
        for vari_index=1:vari_num
            cov=cov+X_dis_sq(:,:,vari_index)*theta(vari_index);
        end
        cov=exp(-cov)+eye(x_num)*1e-6;
        
        % F is base funcion fval of data_point_x
        F_reg=[ones(x_num,1),X];
        
%         % calculation negative log likelihood
%         L=chol(cov)';
%         fval_reg_L=L\fval_reg__;
%         beta__=(fval_reg_L'*fval_reg_L)\(fval_reg_L'*(L\Y));
%         Y_fbeta=Y-fval_reg__*beta__;
%         Y_fbeta_L=L\Y_fbeta;
% %         fval=(Y_fbeta_L'*Y_fbeta_L)/x_num;
%         fval=x_num/2*log((Y_fbeta_L'*Y_fbeta_L)/x_num)+0.5*sum(log(diag(L)));

        % calculation negative log likelihood

        inv_cov=cov\eye(x_num);
        inv_FTRF=(F_reg'*inv_cov*F_reg)\eye(1+vari_num);
        miu=inv_FTRF*(F_reg'*inv_cov*Y);
        Y_Fmiu=Y-F_reg*miu;
        sigma2=(Y_Fmiu'*inv_cov*Y_Fmiu)/x_num;
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
    function [cov,inv_cov,fval_reg,beta,sigma_sq]=interpKriging...
            (X_dis_sq,X,Y,x_number,variable_number,theta)
        % kriging interpolation kernel function
        %
        cov=zeros(x_number,x_number);
        for vari_index=1:variable_number
            cov=cov+X_dis_sq(:,:,vari_index)*theta(vari_index);
        end
        cov=exp(-cov)+eye(x_number)*1e-6;
        
        % F is base funcion fval of data_point_x
%         fval_reg=ones(x_number,1); % zero
        fval_reg=[ones(x_number,1),X]; % linear
        
        % coefficient calculation
        inv_cov=cov\eye(x_number);
        beta=(fval_reg'*inv_cov*fval_reg)\fval_reg'*inv_cov*Y;
        sigma_sq=(Y-fval_reg*beta)'*inv_cov*(Y-fval_reg*beta)/x_number;
    end
    function [Y_pred,Var_pred]=interpKrigingPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            x_num,vari_num,theta,beta,gama,sigma_sq,...
            inv_cov,fval_reg,FTRF,X_pred)
        % kriging interpolation predict function
        % input predict_x and kriging model
        % predict_x is row vector
        % output the predict value
        %
        [x_pred_num,~]=size(X_pred);

        % normalize data
        X_pred=(X_pred-aver_X)./stdD_X;
        
        % predict covariance
        predict_cov=zeros(x_num,x_pred_num);
        for vari_index=1:vari_num
            predict_cov=predict_cov+...
                (X_nomlz(:,vari_index)-X_pred(:,vari_index)').^2*theta(vari_index);
        end
        predict_cov=exp(-predict_cov);

%         predict_fval_reg__=1; % zero
        predict_fval_reg__=[ones(x_pred_num,1),X_pred]; % linear
        Y_pred=predict_fval_reg__*beta+predict_cov'*gama;
        
        % predict variance
        u__=fval_reg'*inv_cov*predict_cov;
        Var_pred=sigma_sq*...
            (1+u__'/FTRF*u__+...
            -predict_cov'*inv_cov*predict_cov);
        
        % normalize data
        Y_pred=Y_pred*stdD_Y+aver_Y;
        Var_pred=diag(Var_pred)*stdD_Y*stdD_Y;
    end
end
