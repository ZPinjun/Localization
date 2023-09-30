function Pout = get_outage_probability(sp,cp,threshold,k,method)

%% get parameters
N = size(sp.Ps_local,2);   % # of SAs
N_Sn = sp.NS_dim(1)*sp.NS_dim(2);   % # of antenna elements in each SA
visi_mat = sp.visi_mat;
Kr = cp.Kr;   % K-factor of Ricean fading
ws = cp.ws;
wbm = cp.wb;
sigman = cp.sigma;
P = cp.P;
H_bar = cp.H_bar;
Gmn = cp.Gmn;
m = cp.m_tilde;


%% calculate outage probability
if method == "analytical"
    Pout = ones(N,1);
    H_bar_k = cell2mat(H_bar(k));
    for n = 1:N
        if visi_mat(m,n) == 1
            H_bar_k_mn = H_bar_k((n-1)*N_Sn+(1:N_Sn),:,m);
            wsn = ws(:,n);
            G_mn = Gmn(m,n,k);
            mu = abs(wsn.' * H_bar_k_mn * wbm);
            sigma2 = norm(G_mn*wbm,2)^2 * norm(wsn,2)^2 / (2*(Kr + 1));
            threshold_new = sqrt( threshold*sigman^2/P );
            pd = makedist('Rician','s',mu,'sigma',sqrt(sigma2));
            Pout(n) = cdf(pd,threshold_new);
        end
    end
elseif method == "empirical"
    repeatNum = 300;
    H_bar_k = cell2mat(H_bar(k));
    SNR = zeros(N,repeatNum);
    for i = 1:repeatNum
        H_tilde = gen_channel_Rayleigh(sp,cp);
        H_tilde_k = cell2mat(H_tilde(k));
        for n = 1:N
            wsn = ws(:,n);
            H_tilde_k_mn = H_tilde_k((n-1)*N_Sn+(1:N_Sn),:,m);
            H_bar_k_mn = H_bar_k((n-1)*N_Sn+(1:N_Sn),:,m);
            H_k_mn = H_bar_k_mn + H_tilde_k_mn;
            SNR(n,i) = P * abs(wsn.'*H_k_mn*wbm)^2 / (sigman^2);
        end
    end
    SNR(SNR>=threshold) = threshold + 1;
    SNR(SNR<threshold) = 1;
    SNR(SNR==threshold+1) = 0;
    Pout = (sum(SNR.')/repeatNum).';
end



end

