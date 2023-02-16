function f = myobj5(x)
%-------Greiwank Problem---------%

persistent o M
[ps,D]=size(x);

load rastrigin_func_data
if length(o)>=D
    o=o(1:D);
else
    o=-5+10*rand(1,D);
end
c=2;
if D==2,load rastrigin_M_D2,
elseif D==10,load rastrigin_M_D10,
elseif D==30,load rastrigin_M_D30,
elseif D==50,load rastrigin_M_D50,
else
    M = rot_matrix(D,c);
end

x=x-repmat(o,ps,1);
x=x*M;
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2) - 330;

end


function M=rot_matrix(D,c)
A=normrnd(0,1,D,D);
P=cGram_Schmidt(A);
A=normrnd(0,1,D,D);
Q=cGram_Schmidt(A);
u=rand(1,D);
D=c.^((u-min(u))./(max(u)-min(u)));
D=diag(D);
M=P*D*Q;
end