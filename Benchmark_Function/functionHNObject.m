function har=functionHNObject(x)
% Hartman function'
%
% variable_number=6;
% object_function=@(x) functionHNObject(x);
% A=[];
% B=[];
% Aeq=[];
% Beq=[];
% low_bou=zeros(variable_number,1);
% up_bou=ones(variable_number,1);
% nonlcon_function=[];
% cheapcon_function=[];
% model_function=[];
% 
% x_best=[0.2017;0.1500;0.4769;0.2753;0.3117;0.6573] fval_best=-3.3224
% 
[row,col]=size(x);
if row~=1&&col~=1
    error('wrong x0 format input, it should be a VECTOR')
end
x=x(:)';
% get coefficient
coe=[1	10	3	17	3.5	1.7	8	1	0.1312	0.1696	0.5569	0.0124	0.8283	0.5886;
2	0.05	10	17	0.1	8	14	1.2	0.2329	0.4135	0.8307	0.3736	0.1004	0.9991;
3	3	3.5	1.7	10	17	8	3	0.2348	0.1451	0.3522	0.2883	0.3047	0.6650;
4	17	8	0.05	10	0.1	14	3.2	0.4047	0.8828	0.8732	0.5743	0.1091	0.0381;];

alpha=coe(:,2:7);
c=coe(:,8);
p=coe(:,9:14);

har=0;
for i=1:4
    hari=0;
    for j=1:6
        hari=alpha(i,j).*(x(:,j)-p(i,j)).^2+hari;
    end
    har=c(i)*exp(-hari)+har;
end
har=-har;

end