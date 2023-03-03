clc;
clear;
close all hidden;

benchmark=BenchmarkFunction();

% % calculate single loss
% variable_number=2;
% sample_number=100;
% x_list=lhsdesign...
%     (sample_number,variable_number);
% fval_list=zeros(sample_number,1);
% for x_index=1:sample_number
%     fval=benchmark.singleBRObject(x_list(x_index,:));
%     fval_list(x_index)=fval(1);
% end
% 
% K=10;
% verify_type='RMSE';
% interp_pre_model_function=@(X,Y) interpKrigingPreModel(X,Y);
% loss=interpCrossVerify...
%     (interp_pre_model_function,K,x_list,fval_list,verify_type)

% compare loss each interp model
variable_number=6;
sample_number=100;
K=10;
verify_type='RMSE';

repeat_time=10;
NRMSE_K_data=zeros(repeat_time,1);
NRMSE_KS_data=zeros(repeat_time,1);

for repeat_index=1:repeat_time
    x_list=lhsdesign...
        (sample_number,variable_number);
    fval_list=zeros(sample_number,1);
    for x_index=1:sample_number
       fval_list(x_index)=benchmark.singleHNObject(x_list(x_index,:));
    end
    
    interp_pre_model_function=@(X,Y) interpKrigingPreModel(X,Y);
    NRMSE_K_data(repeat_index)=interpCrossVerify...
        (interp_pre_model_function,K,x_list,fval_list,verify_type);
    
    interp_pre_model_function=@(X,Y) interpKrigingPreModelM(X,Y);
    NRMSE_KS_data(repeat_index)=interpCrossVerify...
        (interp_pre_model_function,K,x_list,fval_list,verify_type);
end

boxplot([NRMSE_K_data,NRMSE_KS_data]);

