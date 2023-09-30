% meas: angle and delay measurements
% np: known parameters

function Ru_est = estimate_Ru_LS(np)

meas = np.Meas;
Rbs = np.Rbs;
Rsa = np.Rsa;

D = size(meas, 2);
A = zeros(3,D);
B = zeros(3,D);
for i = 1:D
    A(:,i) = -Rbs(:,:,meas(1,i))*[cos(meas(2+2,i)) * cos(meas(1+2,i)); cos(meas(2+2,i)) * sin(meas(1+2,i)); sin(meas(2+2,i))];
    B(:,i) = Rsa(:,:,meas(2,i))*[cos(meas(4+2,i)) * cos(meas(3+2,i)); cos(meas(4+2,i)) * sin(meas(3+2,i)); sin(meas(4+2,i))];
end

ABt = A*B.';
%% =================== use manopt toolbox
% % Create the problem structure.
% manifold = rotationsfactory(3, 1);
% problem.M = manifold;
% % Define the problem cost function and its Euclidean derivatives.
% problem.cost = @(X) -X(:).'*ABt(:);
% problem.egrad = @(X) -ABt;
% problem.ehess = @(X, S) manifold.zerovec(X);
% % Solve the problem
% Ru_est = trustregions(problem);

%% =================== use closed-form solution
% Compare with the known optimal solution
[U, ~, V] = svd(ABt);
UVt = U*V.';
% The determinant of UVt is either 1 or -1, in theory
if abs(1 - det(UVt)) < 1e-10
    Xopt = UVt;
elseif abs(-1 - det(UVt)) < 1e-10
    % UVt is in O(n) but not SO(n). This is easily corrected for:
    J = diag([ones(3-1,1);-1]);
    Xopt = U*J*V.';
else
    error('Should never ever happen ...');
end
Ru_est = Xopt;
    
end

