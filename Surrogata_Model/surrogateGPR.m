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
% [predict_function,kriging_model]=interpGaussPreModel(X,Y);
% check_error_list=zeros(size(matPntHF,1)/K_fold,1);
% for check_index=1:(size(matPntHF,1)/K_fold)
%     x_index=HF_check(check_index);
%     check_error_list(check_index)=(predict_function(matPntHF(x_index,:))-matObjHF(x_index,1));
% end
% disp(['max error: ',num2str(max(check_error_list))]);

% load('PK.mat')

load('matlab.mat');
Y=f;

[predict_function,kriging_model]=interpGaussPreModel(X,Y);
figure_handle=figure(1);
interpVisualize(kriging_model,low_bou,up_bou,figure_handle)

function [predict_function,GPR_model]=interpGaussPreModel...
    (X,Y,theta)
% generate gauss process regression model, version 0
% X, Y are x_number x variable_number matrix
% aver_X,stdD_X is 1 x x_number matrix
% theta beta gama sigma_sq is normalizede, so predict y is normalize
%
% input:
% X, Y (initial data, which are real data)
% theta (hyperparameter, len, eta)
%
% output:
% predict_function, GPR_model (a gauss process regression model)
%
% Copyright 2023.2 Adel
%
[x_number,variable_number]=size(X);
if nargin < 3
    theta=[ones(1,variable_number)*0.5,10];
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
for rank=1:x_number
    for colume=1:rank-1
        X_dis_sq(rank,colume,:)=X_dis_sq(colume,rank,:);
    end
    for colume=rank:x_number
        X_dis_sq(rank,colume,:)=(X_nomlz(rank,:)-X_nomlz(colume,:)).^2;
    end
end

log_likelihood_function=@(theta) logLikelihoodFunction...
    (X_dis_sq,Y,x_number,variable_number,theta);
% optimal to get hyperparameter
object_function=@(x) objectFunctionGauss(log_likelihood_function,x);
low_bou_theta=1e-6*ones(1,variable_number+1);
up_bou_theta=20*ones(1,variable_number+1);

fmincon_option=optimoptions(@fmincon,'Display','iter-detailed',...
    'OptimalityTolerance',1e-2,...
    'FiniteDifferenceStepSize',1e-5,...,
    'MaxIterations',10,'SpecifyObjectiveGradient',true);
[theta,~,~,~]=fmincon...
    (object_function,theta,[],[],[],[],low_bou_theta,up_bou_theta,[],fmincon_option)

len=theta(1:variable_number);
eta=theta(variable_number+1);
% obtain other paramter
[KX,inv_KX,...
    ~]=getCovariance(len,eta,...
    X_dis_sq,x_number,variable_number);

% initialization predict function
predict_function=@(predict_x) interpGaussPredictor...
    (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
    Y_nomlz,len,eta,...
    KX,inv_KX,predict_x);

GPR_model.X=X;
GPR_model.Y=Y;
GPR_model.X_normalize=X_nomlz;
GPR_model.Y_normalize=Y_nomlz;
GPR_model.covariance=KX;
GPR_model.inv_covariance=inv_KX;

GPR_model.theta=theta;
GPR_model.aver_X=aver_X;
GPR_model.stdD_X=stdD_X;
GPR_model.aver_Y=aver_Y;
GPR_model.stdD_Y=stdD_Y;

GPR_model.predict_function=predict_function;

    function [predict_fval,predict_variance]=interpGaussPredictor...
            (X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
            Y,len,eta,...
            KX,inv_KX,predict_x)
        % mu_pre range inf to -inf
        % class is 1 or 0
        % x input is 1 x variable vector
        %
        predict_x=predict_x(:)';
        % normalize data
        predict_x=(predict_x-aver_X)./stdD_X;
        % cov of x and X
        x_cov__=eta*exp(-sum((predict_x-X_nomlz).^2 ./ (len).^2/2,2));

        % get miu and variance of predict x
        predict_fval=x_cov__'*inv_KX*Y;
        predict_variance=eta-x_cov__'*((KX)\x_cov__);

        % normalize data
        predict_fval=predict_fval*stdD_Y+aver_Y;
        predict_variance=predict_variance*stdD_Y*stdD_Y;
    end

    function [fval,gradient,hessian]=logLikelihoodFunction...
            (X_dis_sq,Y,x_number,variable_number,theta)
        % Y is
        % value, X is data
        % hyperparameter is x_number v, variable_number len, eta
        %
        len__=theta(1:variable_number);
        eta__=theta(variable_number+1);

        % obtian covariance matrix
        [K,inv_K,...
            exp_dis__,eta_exp_dis__]=getCovariance...
            (len__,eta__,...
            X_dis_sq,x_number,variable_number);

        detK=det(K);

        fval=-0.5*Y'*inv_K*Y-0.5*log(detK)-x_number/2*log(2*pi);

        if nargout > 1
            % get gradient
            gradient=zeros(variable_number+1,1);

            % var: len1, len2, ... eta
            dK_dvar=zeros(x_number,x_number,variable_number+1);
            dinv_K_dvar=zeros(x_number,x_number,variable_number+1);
            ddetK_dvar=zeros(1,1,variable_number+1);
            % len
            for len_index__=1:variable_number
                dK_dlen=eta_exp_dis__.*X_dis_sq(:,:,len_index__)/len__(len_index__)^3;
                dinv_K_dlen=-inv_K*dK_dlen*inv_K;
                ddetK_dlen=detK*trace(inv_K*dK_dlen);

                dK_dvar(:,:,len_index__)=dK_dlen;
                dinv_K_dvar(:,:,len_index__)=dinv_K_dlen;
                ddetK_dvar(:,:,len_index__)=ddetK_dlen;

                gradient(len_index__)=-0.5*Y'*dinv_K_dlen*Y+...
                    -0.5/detK*ddetK_dlen;
            end
            % eta
            dK_deta=exp_dis__;
            dinv_K_deta=-inv_K*dK_deta*inv_K;
            ddetK_deta=detK*trace(inv_K*dK_deta);

            dK_dvar(:,:,end)=dK_deta;
            dinv_K_dvar(:,:,end)=dinv_K_deta;
            ddetK_dvar(:,:,end)=ddetK_deta;

            gradient(variable_number+1)=-0.5*Y'*dinv_K_deta*Y+...
                -0.5/detK*ddetK_deta;
        end

        if nargout > 2
            % get hessian
            % var: len1, len2, ... eta
            hessian=zeros(variable_number+1);
            for len_i=1:variable_number
                for len_j=1:len_i-1
                    hessian(len_i,len_j)=hessian(len_j,len_i);
                end
                len_j=len_i;
                ddK_dii=dK_dvar(:,:,len_i).*X_dis_sq(:,:,len_j)/len__(len_j)^3+...
                    eta_exp_dis__.*(-3*X_dis_sq(:,:,len_j)/len__(len_j)^4);
                ddinv_K_dij=-dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)*inv_K+...
                    -inv_K*ddK_dii*inv_K+...
                    -inv_K*dK_dvar(:,:,len_i)*dinv_K_dvar(:,:,len_j);
                dddetK_dij=ddetK_dvar(:,:,len_j)*trace(inv_K*dK_dvar(:,:,len_i))+...
                    detK*trace(dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)+inv_K*ddK_dii);
                hessian(len_i,len_j)=-0.5*Y'*ddinv_K_dij*Y-0.5/detK*dddetK_dij+...
                    0.5/detK^2*ddetK_dvar(:,:,len_j)*ddetK_dvar(:,:,len_i);

                for len_j=len_i+1:variable_number
                    ddK_dij=dK_dvar(:,:,len_i).*X_dis_sq(:,:,len_j)/len__(len_j)^3;
                    ddinv_K_dij=-dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)*inv_K+...
                        -inv_K*ddK_dij*inv_K+...
                        -inv_K*dK_dvar(:,:,len_i)*dinv_K_dvar(:,:,len_j);
                    dddetK_dij=ddetK_dvar(:,:,len_j)*trace(inv_K*dK_dvar(:,:,len_i))+...
                        detK*trace(dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)+inv_K*ddK_dij);
                    hessian(len_i,len_j)=-0.5*Y'*ddinv_K_dij*Y-0.5/detK*dddetK_dij+...
                        0.5/detK^2*ddetK_dvar(:,:,len_j)*ddetK_dvar(:,:,len_i);
                end

                len_j=variable_number+1; % eta
                ddK_dieta=exp_dis__.*X_dis_sq(:,:,len_i)/len__(len_i)^3;
                ddinv_K_dij=-dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)*inv_K+...
                    -inv_K*ddK_dieta*inv_K+...
                    -inv_K*dK_dvar(:,:,len_i)*dinv_K_dvar(:,:,len_j);
                dddetK_dij=ddetK_dvar(:,:,len_j)*trace(inv_K*dK_dvar(:,:,len_i))+...
                    detK*trace(dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)+inv_K*ddK_dieta);
                hessian(len_i,len_j)=-0.5*Y'*ddinv_K_dij*Y-0.5/detK*dddetK_dij+...
                    0.5/detK^2*ddetK_dvar(:,:,len_j)*ddetK_dvar(:,:,len_i);
            end

            % eta eta
            len_i=variable_number+1;
            for len_j=1:len_i-1
                hessian(len_i,len_j)=hessian(len_j,len_i);
            end
            len_j=len_i;
            ddK_detaeta=0;
            ddinv_K_dij=-dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)*inv_K+...
                -inv_K*ddK_detaeta*inv_K+...
                -inv_K*dK_dvar(:,:,len_i)*dinv_K_dvar(:,:,len_j);
            dddetK_dij=ddetK_dvar(:,:,len_j)*trace(inv_K*dK_dvar(:,:,len_i))+...
                detK*trace(dinv_K_dvar(:,:,len_j)*dK_dvar(:,:,len_i)+inv_K*ddK_detaeta);
            hessian(len_i,len_j)=-0.5*Y'*ddinv_K_dij*Y-0.5/detK*dddetK_dij+...
                0.5/detK^2*ddetK_dvar(:,:,len_j)*ddetK_dvar(:,:,len_i);

        end
    end

    function [fval,gradient]=objectFunctionGauss(initial_function,x)
        [fval,gradient]=initial_function(x);
        fval=-fval;
        gradient=-gradient';
    end

    function [KX__,inv_KX__,...
            exp_dis__,eta_exp_dis__]=getCovariance(len,eta,...
            X_dis_sq,x_number,variable_number)
        % obtain covariance of x
        %
        % exp of x__x with theta
        exp_dis__=zeros(x_number);
        for rank_index__=1:x_number
            % symmetry
            for colume_index__=1:rank_index__-1
                exp_dis__(rank_index__,colume_index__)=exp_dis__(colume_index__,rank_index__);
            end

            % diagonal
            exp_dis__(rank_index__,rank_index__)=1+1e-6;

            % initial
            for colume_index__=rank_index__+1:x_number
                temp__=X_dis_sq(rank_index__,colume_index__,:);
                exp_dis__(rank_index__,colume_index__)=...
                    exp(-(0.5./len.^2)*temp__(:));
            end
        end

        KX__=eta*exp_dis__;
        inv_KX__=inv(KX__);

        eta_exp_dis__=KX__;
    end
end
