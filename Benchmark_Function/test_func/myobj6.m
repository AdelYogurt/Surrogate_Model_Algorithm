function fit = myobj6(x)

persistent  fun_num func o sigma lamda bias M

[ps,D]=size(x);
initial_flag=1;
fun_num=10;
load hybrid_func2_data % saved the predefined optima
if length(o(1,:))>=D
    o=o(:,1:D);
else
    o=-5+10*rand(fun_num,D);
end
o(10,:)=0;
func.f1=str2func('fackley');
func.f2=str2func('fackley');
func.f3=str2func('frastrigin');
func.f4=str2func('frastrigin');
func.f5=str2func('fsphere');
func.f6=str2func('fsphere');
func.f7=str2func('fweierstrass');
func.f8=str2func('fweierstrass');
func.f9=str2func('fgriewank');
func.f10=str2func('fgriewank');
bias=((1:fun_num)-1).*100;
sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
lamda=repmat(lamda,1,D);
c=[2 3 2 3 2 3 20 30 200 300];
if D==2,load hybrid_func2_M_D2,
elseif D==10,load hybrid_func2_M_D10,
elseif D==30,load hybrid_func2_M_D30,
elseif D==50,load hybrid_func2_M_D50,
else
    for i=1:fun_num
        eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
    end
end

fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
end