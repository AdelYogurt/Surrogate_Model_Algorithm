function problem_info = problem_define(dimension, problem)

%The test problems are showed as following
% Ellipsoid; Rosenbrock; Ackely; Greiwank; 
% Shifted Rotated Rastrigin;
% Rotated hybrid Composition Function

switch problem
    
    case 'Ellipsoid'
        objfun = @myobj1;
        lb = -5.12 * ones(1, dimension);
        ub = 5.12 * ones(1, dimension);
        
    case 'Rosenbrock'
        objfun = @myobj2;
        lb = -2.048 * ones(1, dimension);
        ub = 2.048 * ones(1, dimension);
        
    case 'Ackely'
        objfun = @myobj3;
        lb = -32.768 * ones(1, dimension);
        ub = 32.768 * ones(1, dimension);
        
    case 'Greiwank'
        objfun = @myobj4;
        lb = -600 * ones(1, dimension);
        ub = 600 * ones(1, dimension);
        
    case 'Shifted Rotated Rastrigin'
        objfun = @myobj5;
        lb = -5 * ones(1, dimension);
        ub = 5 * ones(1, dimension);
        
    case 'Rotated Hybrid Composition Function'
        objfun = @myobj6;
        lb = -5 * ones(1, dimension);
        ub = 5 * ones(1, dimension);
end

problem_info.objfun = objfun;
problem_info.lb = lb;
problem_info.ub = ub;
problem_info.dimension = dimension;

end