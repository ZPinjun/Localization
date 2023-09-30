function [RMSE_Pu_LS,RMSE_Ru_LS,RMSE_Pu_ML,RMSE_Ru_ML] = get_RMSE(repeatNum,cp,sp,Pu,Ru,CCRB)

% get parameters
visi_mat = sp.visi_mat;
D = sum(sum(visi_mat));

% get noise covariance matrix
FIM_eta_ch = get_FIM_channel_params(cp,sp);   
A = FIM_eta_ch(1:5*D,1:5*D);
B = FIM_eta_ch(1:5*D, (5*D+1):7*D);
C = FIM_eta_ch((5*D+1):7*D,(5*D+1):7*D);
Sigma = (A - B*C^(-1)*B.')^(-1); % measurement covariance matrix

% get measurements
np.Meas = sp.Meas;
np.Meas(:,isnan(np.Meas(3,:))) = [];

% % Add noise
meas1 = np.Meas(3:end,:);
meas1(1:(end-1),:) = deg2rad(meas1(1:end-1,:));
mu = reshape(meas1.',1,[]);
measn_all = mvnrnd(mu,Sigma,repeatNum);

np.Rbs = sp.Rb;    % Ori of BSs (global coordinate system)
np.Rsa = sp.Rs_local;   % Ori of SAs (UE's local coordinate system)
np.pbs = sp.Pb;     % Pos of BSs (global coordinate system)
np.dsa = sp.Ps_local; % Pos of SAs (UE's local coordinate system)
np.Sigma = Sigma;

data_LS = zeros(2,repeatNum);
data_ML = zeros(2,repeatNum);
for i = 1:repeatNum
    %% call algorithm
    % ========================== add noise
    measn = measn_all(i,:);
    measn = (reshape(measn, [], 5)).';
    np.Meas(3:end,:) = measn;
    % ========================== call algorithm  
    % LS
    Ru_est_LS = estimate_Ru_LS(np);
    r14_est_LS = estimate_pu_rho_LS(Ru_est_LS, np);
    X_init.r14 = r14_est_LS;
    X_init.Ru = Ru_est_LS;
    % ML
    [Pu_est, ~, Ru_est] = localization_algorithm(np,X_init);
    % save data
    data_LS(1,i) = norm(r14_est_LS(1:3) - Pu, 2)^2;
    data_LS(2,i) = norm(Ru_est_LS(:) - Ru(:), 2)^2;
    data_ML(1,i) = norm(Pu_est - Pu, 2)^2;
    data_ML(2,i) = norm(Ru_est(:) - Ru(:), 2)^2;

    % show results
    figure(10)
    semilogy(1:repeatNum,ones(1,repeatNum)*CCRB(1)^2,'k-','LineWidth',2); hold on 
    semilogy(1:repeatNum,ones(1,repeatNum)*(CCRB(2)/40.5142)^2,'b-','LineWidth',2); hold on 


    semilogy(1:repeatNum,data_LS(1,:),'r^-'); hold on
    semilogy(1:repeatNum,data_LS(2,:),'g^-'); hold on
    semilogy(1:repeatNum,data_ML(1,:),'rv-'); hold on
    semilogy(1:repeatNum,data_ML(2,:),'gv-'); hold on

    semilogy(1:repeatNum,ones(1,repeatNum)*mean(data_LS(1,1:i)),'r:'); hold on
    semilogy(1:repeatNum,ones(1,repeatNum)*mean(data_LS(2,1:i)),'g:'); hold on
    semilogy(1:repeatNum,ones(1,repeatNum)*mean(data_ML(1,1:i)),'r--'); hold on
    semilogy(1:repeatNum,ones(1,repeatNum)*mean(data_ML(2,1:i)),'g--'); hold off

end

RMSE_Pu_LS = sqrt(mean(data_LS(1,:).'));
RMSE_Ru_LS = 40.5142*sqrt(mean(data_LS(2,:).'));
RMSE_Pu_ML = sqrt(mean(data_ML(1,:).'));
RMSE_Ru_ML = 40.5142*sqrt(mean(data_ML(2,:).'));

end

