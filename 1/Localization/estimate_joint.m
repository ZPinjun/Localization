function X_est = estimate_joint(X_init, np)

%% get parameters
meas = np.Meas;
Sigma = np.Sigma;
Rbs = np.Rbs;
Rsa = np.Rsa;
pbs = np.pbs;
dsa = np.dsa;



%% =================== use manopt toolbox
% Create the problem structure
so3Plus_manifold = productmanifold(struct( ...
    'r14', euclideanfactory(4,1), ...
    'Ru', rotationsfactory(3) ...
));
problem.M = so3Plus_manifold;
% Define the problem cost function and its Euclidean derivatives
problem.cost = @(X) cost_function(X, meas, Sigma, Rbs, Rsa, pbs, dsa);
problem.grad = @(X) (so3Plus_manifold.egrad2rgrad(X, egrad_function(X, meas, Sigma, Rbs, Rsa, pbs, dsa)));
% Solve the problem
opt.minstepsize = 1e-20;
opt.maxiter = 5000;
% opt.verbosity = 0;
[X_est,~] = steepestdescent(problem, X_init, opt);







function f_cost = cost_function(X, meas, Sigma, Rbs, Rsa, pb, d_n)
    pu = X.r14(1:3);
    rho = X.r14(4);
    Ru = X.Ru;
    eta_hat = reshape(meas(3:end,:).',[],1);
    % eta
    eta = meas;
    for i = 1:size(eta,2)
        m = meas(1,i);   % label of BS
        n = meas(2,i);   % label of SA
        PBS = pb(:,m);
        PUE = pu + Ru*d_n(:,n);
        % AOD
        tBU = PUE - PBS;   % in the global coordinate system
        tBU_BSCS = Rbs(:,:,m).'*tBU;   % in the BS's local coordinate system  
        AOD = get_angle_from_dir(tBU_BSCS);
        % AOA
        tUB = PBS - PUE;   % in the global coordinate system
        tUB_UECS = (Ru*Rsa(:,:,n)).'*tUB;   % in the SA's local coordinate system  
        AOA = get_angle_from_dir(tUB_UECS);
        % delay, in [mm], note rho is in [m]
        tau = norm(tUB, 2) + rho;
        % get eta
        eta(3:end,i) = [deg2rad(AOD); deg2rad(AOA); tau];
    end
    eta = reshape(eta(3:end,:).',[],1);
    f_cost = 0.5*(eta_hat - eta).'*Sigma^(-1)*(eta_hat - eta);
end

function egrad = egrad_function(X, meas, Sigma, Rbs, Rsa, pb, d_n)
    pu = X.r14(1:3);
    rho = X.r14(4);
    Ru = X.Ru;
    eta_hat = reshape(meas(3:end,:).',[],1);
    D = size(meas,2);
    % eta
    eta = meas;
    for i = 1:D
        m = meas(1,i);   % label of BS
        n = meas(2,i);   % label of SA
        PBS = pb(:,m);
        PUE = pu + Ru*d_n(:,n);
        % AOD
        tBU = PUE - PBS;   % in the global coordinate system
        tBU_BSCS = Rbs(:,:,m).'*tBU;   % in the BS's local coordinate system  
        AOD = get_angle_from_dir(tBU_BSCS);
        % AOA
        tUB = PBS - PUE;   % in the global coordinate system
        tUB_UECS = (Ru*Rsa(:,:,n)).'*tUB;   % in the SA's local coordinate system  
        AOA = get_angle_from_dir(tUB_UECS);
        % delay, in [mm]
        tau = norm(tUB, 2) + rho;
        % get eta
        eta(3:end,i) = [deg2rad(AOD); deg2rad(AOA); tau];
    end
    eta = reshape(eta(3:end,:).',[],1);
    % df_deta
    df_deta = -Sigma^(-1)*(eta_hat - eta);
    % deta_dr14
    deta_dr14 = zeros(5*D,4);
    deta_dRu = zeros(5*D,9);
    u1 = [1 0 0].';
    u2 = [0 1 0].';
    u3 = [0 0 1].';
    for j = 1:D
        m = meas(1,j);   % label of BS
        n = meas(2,j);   % label of SA
        % SA
        dn = d_n(:,n);
        Rn = Rsa(:,:,n);
        % BS
        pBm = pb(:,m);
        RBm = Rbs(:,:,m);
        %
        tbs = pu + Ru*dn - pBm;
        normtbs = norm(tbs,2);
        dir_BS = RBm.'*tbs;
        dir_SA = Rn.'*Ru.'*tbs;
        termdn = [dn(1)*eye(3); dn(2)*eye(3); dn(3)*eye(3)];
        % row1
        dAODazi_dpu = (dir_BS(1)*RBm*u2 - dir_BS(2)*RBm*u1) / (dir_BS(1)^2+dir_BS(2)^2);
        dAODazi_dRu = (dir_BS(1)*RBm*u2*dn.' - dir_BS(2)*RBm*u1*dn.') / (dir_BS(1)^2+dir_BS(2)^2);
        dAODazi_dRu_vec = reshape(dAODazi_dRu, [9,1]);
        % row2
        term2 = 1/sqrt(1 - (dir_BS(3)/normtbs)^2);
        dAODele_dpu = term2 * (  RBm*u3/normtbs - dir_BS(3)*tbs/normtbs^3  );
        dAODele_dRu_vec = term2 * (  reshape(RBm*u3*dn.', [9,1])/normtbs - dir_BS(3)*termdn*tbs/normtbs^3  );
        % row3
        dAOAazi_dpu = (dir_SA(1)*Ru*Rn*u2 - dir_SA(2)*Ru*Rn*u1)/(dir_SA(1)^2+dir_SA(2)^2);
        dtermu1 =  (pu - pBm)*u2.'*Rn.' + Ru*(Rn*u2*dn.' + dn*u2.'*Rn.');
        dtermu2 =  (pu - pBm)*u1.'*Rn.' + Ru*(Rn*u1*dn.' + dn*u1.'*Rn.');
        dAOAazi_dRu = (dir_SA(1)*dtermu1 - dir_SA(2)*dtermu2)/(dir_SA(1)^2+dir_SA(2)^2);
        dAOAazi_dRu_vec = reshape(dAOAazi_dRu, [9,1]);
        % row4
        term4 = 1/sqrt(1 - (dir_SA(3)/normtbs )^2);
        dAOAele_dpu = -term4 * ( Ru*Rn*u3/normtbs - dir_SA(3)*tbs/normtbs^3 );
        dtermu =  (pu - pBm)*u3.'*Rn.' + Ru*(Rn*u3*dn.' + dn*u3.'*Rn.');
        dAOAele_dRu_vec = -term4 * (   reshape(dtermu, [9,1])/normtbs - dir_SA(3)*termdn*tbs/normtbs^3   );
        % row5
        dtau_dpu = tbs/normtbs;
        dtau_dRu_vec = termdn*tbs/normtbs;
        
        dl_dpu = [dAODazi_dpu.'; dAODele_dpu.'; dAOAazi_dpu.'; dAOAele_dpu.'; dtau_dpu.'];
        deta_dr14(([0,D,2*D,3*D,4*D]+j),:) = [dl_dpu,[0;0;0;0;1]];
        deta_dRu(([0,D,2*D,3*D,4*D]+j),:) = [dAODazi_dRu_vec.'; dAODele_dRu_vec.'; dAOAazi_dRu_vec.'; dAOAele_dRu_vec.'; dtau_dRu_vec.'];
    end
    egrad.r14 = deta_dr14.'*df_deta;
    egrad.Ru = reshape(deta_dRu.'*df_deta,[3,3]);
end

end

