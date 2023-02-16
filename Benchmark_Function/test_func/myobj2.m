function f = myobj2(x)
%-------Rosenbrock Problem---------%

xi = x;
xi(size(x,2)) = [];
xi_1 = x;
xi_1(1) = [];

f = sum(100*(xi.^2 - xi_1).^2 + ( xi - 1 ).^2 );

end
