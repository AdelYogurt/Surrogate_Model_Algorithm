function fval=functionG01ObjectLow(x)
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
fval=functionG01Object(x)*0.9+0.5;
end
