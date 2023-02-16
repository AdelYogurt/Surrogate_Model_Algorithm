function f = myobj4(x)
%-------Shifted Rotated Rastrigin Problem---------%

i = 1 : size(x,2);

f = sum(x.^2/4000) - prod( cos(x/sqrt(i)) ) + 1;

end
