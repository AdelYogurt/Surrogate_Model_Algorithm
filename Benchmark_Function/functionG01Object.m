function f=functionG01Object(x)
% object_function=@functionG01Object;
% object_function_low=@functionG01ObjectLow;
% A=[ 2   2   0   0   0   0   0   0   0   1   1   0   0;
%     2   0   2   0   0   0   0   0   0   1   0   1  0;
%     0   2   2   0   0   0   0   0   0   0   1   1  0;
%     -8  0   0   0   0   0   0   0   0   1   0   0   0;
%     0   -8  0   0   0   0   0   0   0   0   1   0   0;
%     0   0   0   -2  -1  0   0   0   0   1   0   0   0;
%     0   0   0   0   0   -2  -1  0   0   0   1   0   0;
%     0   0   0   0   0   0   0   -2  -1  0   0   1   0;
%     ];
% B=[10;10;10;0;0;0;0;0];
% Aeq=[];
% Beq=[];
% low_bou=zeros(13,1);
% up_bou=ones(13,1);
% up_bou(10:12)=100;
% nonlcon_function=[];
%
sigma1=0;
for i=1:4
    sigma1=sigma1+x(i);
end
sigma2=0;
for i=1:4
    sigma2=sigma2+x(i)^2;
end
sigma3=0;
for i=5:13
    sigma3=x(i)+sigma3;
end
f=5*sigma1-5*sigma2-sigma3;
end
