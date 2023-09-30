function FIM = get_FIM_channel_params(cp, sp)

%% get parameters
N_Bm = sp.NB_dim(1)*sp.NB_dim(2);   % # of antenna elements in each BS
N_Sn = sp.NS_dim(1)*sp.NS_dim(2);   % # of antenna elements in each SA
visi_mat = sp.visi_mat;
Meas = sp.Meas;
Pb = sp.Pb;
Rb = sp.Rb;
Ps = sp.Ps;
Rs = sp.Rs;
pAE_B = sp.pAE_B;
pAE_S = sp.pAE_S;
K = cp.K;   % # of subcarriers
G = cp.G;   % # of transmissions
P = cp.P;
Kr = cp.Kr;
sigma = cp.sigma;
x = cp.x;
wb = cp.wb;
ws = cp.ws;
fk = cp.fk;
fc = cp.fc;
c = cp.c;
Gmn = cp.Gmn;

D = sum(sum(visi_mat));
Meas(:,logical(reshape(1-visi_mat.',1,[]))) = [];



%% get phases of the beamsteering vector in [m]
PsiB = zeros(N_Bm,D);
PsiS = zeros(N_Sn,D);
for i = 1:D
    m = Meas(1,i);
    n = Meas(2,i);
    % AOD
    tAOD_global = Ps(:,n) - Pb(:,m);   % in the global coordinate system
    tAOD_local = Rb(:,:,m).'*tAOD_global;   % in the BS's local coordinate system
    %tAOD_local = normalize(tAOD_local, 'norm', 2);
    tAOD_local = tAOD_local/norm(tAOD_local,2);
    PsiB(:,i) = pAE_B.' * tAOD_local;
    % AOA
    tAOA_global = Pb(:,m) - Ps(:,n);   % in the global coordinate system
    tAOA_local = Rs(:,:,n).'*tAOA_global;   % in the BS's local coordinate system
    %tAOA_local = normalize(tAOA_local, 'norm', 2);
    tAOA_local = tAOA_local/norm(tAOA_local,2);
    PsiS(:,i) = pAE_S.' * tAOA_local;
end



%% calculate FIM
FIM = zeros(7*D, 7*D);
% loop over transmissions
for g = 1:G
    % loop over subcarries
    for k = 1:K
        Sigma_tilde = zeros(D,D);
        D_m_eta_ch = zeros(D,7*D);
        for i = 1:D
            % calculate Sigma_i
            m = Meas(1,i);
            n = Meas(2,i);
            G_mn = Gmn(m, n, k);
            Sigma_tilde(i,i) = sigma^2 + P*G_mn^2/(1+Kr);   % only when ||w_Bm|| = ||w_Sn|| = ||x|| = 1
            % calculate derivative
            x_m = x(m,k,g);
            ws_n = ws(:,n,g);
            wb_m = wb(:,m,g);  
            qmn = sqrt(P)*ws_n.';
            xmn = wb_m*x_m;
            hmn_a = G_mn*sqrt(Kr/(1+Kr));
            aBm = exp(1j*2*pi*fk(k)/c*PsiB(:,i));
            aSn = exp(1j*2*pi*fk(k)/c*PsiS(:,i));
            eta = Meas(3:end,i);
            D_eta_i = get_D_eta(qmn,xmn,hmn_a,aBm,aSn,fk(k),fc,c,eta,pAE_B,pAE_S);
            D_m_eta_ch(i,([0,D,2*D,3*D,4*D,5*D,6*D]+i)) = D_eta_i;
        end
        FIM_gk = D_m_eta_ch' * Sigma_tilde^(-1) * D_m_eta_ch;
        FIM = FIM + 2*real(FIM_gk);
    end
end


end

