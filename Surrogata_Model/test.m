clc;
clear;
close all hidden;

load('matlab.mat');

step=1e-5;

cov__=zeros(x_number,x_number);
for rank_index__=1:x_number
    for colume_index__=1:rank_index__-1
        cov__(rank_index__,colume_index__)=...
            cov__(colume_index__,rank_index__);
    end
    cov__(rank_index__,rank_index__)=1; % stabilize
    for colume_index__=rank_index__+1:x_number
        temp__=X_dis_sq(rank_index__,colume_index__,:);
        cov__(rank_index__,colume_index__)=...
            exp(-temp__(:)'*theta(:));
    end
end

cov_plus__=zeros(x_number,x_number);
for rank_index__=1:x_number
    for colume_index__=1:rank_index__-1
        cov_plus__(rank_index__,colume_index__)=...
            cov_plus__(colume_index__,rank_index__);
    end
    cov_plus__(rank_index__,rank_index__)=1; % stabilize
    for colume_index__=rank_index__+1:x_number
        temp__=X_dis_sq(rank_index__,colume_index__,:);
        cov_plus__(rank_index__,colume_index__)=...
            exp(-temp__(:)'*(theta(:)+[step;0]));
    end
end


inv_fTcovf=1/(fval_reg__'*inv(cov__)*fval_reg__);
beta__=inv_fTcovf*fval_reg__'*inv(cov__)*Y;
Y_fbeta=Y-fval_reg__*beta__;
fval=(Y_fbeta'*inv(cov__)*Y_fbeta)/x_number;

inv_fTcovf_plus=1/(fval_reg__'*inv(cov_plus__)*fval_reg__);
beta_plus__=inv_fTcovf_plus*fval_reg__'*inv(cov_plus__)*Y;
Y_fbeta_plus=Y-fval_reg__*beta_plus__;
fval_plus=(Y_fbeta_plus'*inv(cov_plus__)*Y_fbeta_plus)/x_number;

(inv_fTcovf_plus-inv_fTcovf)/step
(beta_plus__-beta__)/step
% (fval_plus-fval)/step

(inv(cov_plus__)-inv(cov__))/step;
dinv_cov_dtheta=...
    inv_cov__*(X_dis_sq(:,:,variable_index).*cov__)*inv_cov__;

dinv_fTcovf_dtheta=-inv_fTcovf*...
    (fval_reg__'*dinv_cov_dtheta*fval_reg__)*...
    inv_fTcovf
dbeta_dtheta=dinv_fTcovf_dtheta*fval_reg__'*inv_cov__*Y+...
    inv_fTcovf*fval_reg__'*dinv_cov_dtheta*Y


% dy_fbeta_dtheta=Y-fval_reg__*dbeta_dtheta;
% (dy_fbeta_dtheta'*inv_cov__*Y_fbeta+...
%     Y_fbeta'*dinv_cov_dtheta*Y_fbeta+...
%     Y_fbeta'*inv_cov__*dy_fbeta_dtheta)/x_number