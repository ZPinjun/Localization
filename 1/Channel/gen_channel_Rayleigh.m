function H_tilde = gen_channel_Rayleigh(sp, cp)

%% get parameters
M = size(sp.Pb,2);   % # of BSs
N = size(sp.Ps_local,2);   % # of SAs
N_Bm = sp.NB_dim(1)*sp.NB_dim(2);   % # of antenna elements in each BS
N_Sn = sp.NS_dim(1)*sp.NS_dim(2);   % # of antenna elements in each SA
K = cp.K;   % # of subcarriers
Kr = cp.Kr;   % K-factor of Ricean fading
Gmn = cp.Gmn;

%% calculate the scattering component of the channel matrix
H_tilde = cell(1,K);
for k = 1:K
    H_tilde_k = zeros(N*N_Sn, N_Bm, M);
    for m = 1:M
        for n = 1:N
            if Gmn(m,n,k) ~= 0
                H_tilde_k((n-1)*N_Sn+(1:N_Sn),:,m) = Gmn(m,n,k) * randn(N_Sn, N_Bm);
            end
        end
    end
    H_tilde(k) = mat2cell(H_tilde_k/sqrt(1+Kr), N_Sn*N, N_Bm, M);
end

end

