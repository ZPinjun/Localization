function [X_hat,supp] = func_SOMP(A,S,sThres)
% ----------- INPUTS:
% A: sensing matrix
% S: observation matrix
% sThres: sparsity threshold
% ----------- OUTPUTS:
% X_hat: estimated unknown vectors
% supp: the support of the unknown sparse vectors

% normalize the sensing matrix
M = size(A,1);
N = size(A,2);
L = size(S,2);
d = sqrt(sum(A.*conj(A)));
A = A ./ kron( ones(M,1), d );

% initialization
Rk = S;
supp = [];

% iterations
for k = 1: sThres
    [~, indx] = max(sum(abs(A'*Rk),2));
    supp = [supp, indx];
    Ak = A(:,supp);
    Xk = (Ak'*Ak)^(-1)*Ak'*S;
    Bk = Ak*Xk;
    Rk = S - Bk;
end

% form outputs
X_hat = zeros(N,L);
X_hat(supp,:) = Xk;

end

