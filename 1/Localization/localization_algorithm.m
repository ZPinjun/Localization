% Pu: UE position estimate, 3*1, in [m]
% rho: clock offset estmate in [m]
% Ru: UE orientation estimate, 3*3 rotation matrix

function [Pu, rho, Ru] = localization_algorithm(np,X_init)


%% initialization (LS)
if ~exist('X_init','var')
    Ru_est_LS = estimate_Ru_LS(np);
    r14_est_LS = estimate_pu_rho_LS(Ru_est_LS, np);
    X_init.r14 = r14_est_LS;
    X_init.Ru = Ru_est_LS;
end

%% joint manifold optimization
X_est = estimate_joint(X_init, np);
Pu = X_est.r14(1:3);
rho = X_est.r14(4);
Ru = X_est.Ru;
end

