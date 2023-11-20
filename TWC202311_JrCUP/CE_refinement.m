function Xest = CE_refinement(Xinit, y, sp)

%% =================== use manopt toolbox
% Create the problem structure.
manifold = euclideanfactory(8,1);
problem.M = manifold;

% Define the problem cost function and its Euclidean derivatives
problem.cost = @(X) CE_cost_fun(X, y, sp);
%problem.grad = @(X) (manifold.egrad2rgrad(X, egrad_function(X, y, sp)));
% Solve the problem
opt.minstepsize = 1e-10;
opt.maxiter = 40;
opt.tolgradnorm = 1e-9;
opt.verbosity = 0;
warning('off', 'manopt:getGradient:approx');
warning('off', 'manopt:getHessian:approx');
[Xest, ~, ~, ~] = trustregions(problem, Xinit, opt);


end