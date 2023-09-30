function C = get_ergodic_capacity(sp,cp,repeatNum)

%% get parameters
K = cp.K;   % # of subcarriers
P = cp.P;
sigma2 = cp.sigma^2;
visi_mat = sp.visi_mat;
H_bar = cp.H_bar;
wb = cp.wb;   % precoder
ws = cp.ws;   % combiner
N = size(ws,2);   % # of SAs
N_Sn = size(ws,1);   % # of antenna elements in each SA
m_tilde = cp.m_tilde;
df = cp.BW/cp.K;


%% calculate ergodic capacity
C = zeros(K,repeatNum);
m = m_tilde;
if sum(visi_mat(m,:)) > 1e-10   % no need to compute if BS m is invisible
    for i = 1:repeatNum
        H_tilde = gen_channel_Rayleigh(sp,cp);
        for k = 1:K
            % precoder
             H_bar_k = cell2mat(H_bar(k));
             H_tilde_k = cell2mat(H_tilde(k));
            for n = 1:N
                if sum(visi_mat(m,n)) > 1e-10   % no need to compute if SA n is invisible
                    H_bar_m_n = H_bar_k((n-1)*N_Sn+(1:N_Sn),:,m);
                    H_tilde_m_n = H_tilde_k((n-1)*N_Sn+(1:N_Sn),:,m);
                    H_m_n = H_bar_m_n + H_tilde_m_n;
                    SNRm = P*abs(ws(:,n).'*H_m_n*wb)^2/sigma2;
                    % compute sum-rate
                    C(k,i) = C(k,i) + df*log2(1+SNRm);
                end
            end
        end
    end
end

C = sum(sum(C))/repeatNum;


end

