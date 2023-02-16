function f = myobj1(x)
%-------Ellipsoid Problem---------%

dim = size(x, 2);
x2 = x.^2;

martix_i =[1: dim]';

f = x2 * martix_i;

end
