function [x_hat,supp] = func_OMP(A,b,sThres)
% ----------- INPUTS:
% A: sensing matrix
% b: observation vector
% sThres: sparsity threshold
% ----------- OUTPUTS:
% x_hat: estimated unknown vector
% supp: the support of the unknown sparse vector

% normalize the sensing matrix
M = size(A,1);
N = size(A,2);
d = sqrt(sum(A.*conj(A)));
A = A ./ kron( ones(M,1), d );

% initialization
rk = b;
supp = [];

% iterations
for k = 1: sThres
    [~, indx] = max(abs(A'*rk));
    supp = [supp, indx];
    Ak = A(:,supp);
    xk = (Ak'*Ak)^(-1)*Ak'*b;
    bk = Ak*xk;
    rk = b - bk;
end

% form outputs
x_hat = zeros(N,1);
x_hat(supp) = xk;


end

