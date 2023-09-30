function cp = gen_precoder_combiner(sp,cp,pcp)

%% get parameters
K = cp.K;   % # of subcarriers
G = cp.G;   % # of transmissions
M = size(sp.Pb,2);   % # of BSs
N = size(sp.Ps_local,2);   % # of SAs
N_Bm = sp.NB_dim(1)*sp.NB_dim(2);   % # of antenna elements in each BS
N_Sn = sp.NS_dim(1)*sp.NS_dim(2);   % # of antenna elements in each SA
H_bar = cp.H_bar;
H_tilde = cp.H_tilde;
visi_mat = sp.visi_mat;
sigma2 = (cp.sigma)^2;
P = cp.P;
df = cp.BW/cp.K;

if strcmp(pcp.func,'localization')  
    % generate transmitted symbols (before precoder), index: [BS, subcarrier, transmission]
    cp.x = exp(1j*2*pi*rand(M,K,G));
    % generate precoders, index: [AE, BS, transmission]
    cp.wb = exp(1j*2*pi*rand(N_Bm,M,G))/sqrt(N_Bm);
    % generate combiners, index: [AE, SA, transmission]
    cp.ws = exp(1j*2*pi*rand(N_Sn,N,G))/sqrt(N_Sn);
elseif strcmp(pcp.func,'communication')
    % generate transmitted symbols (before precoder), index: [BS, subcarrier, transmission]
    cp.x = exp(1j*2*pi*rand(M,K,G));
    % select a BS with the maximum sum-rate
    data_sumrate = zeros(1,M);
    for m = 1:M
        if sum(visi_mat(m,:)) > 1e-10
            for k = 1:K
                % compute precoder
                H_bar_k = cell2mat(H_bar(k));
                H_tilde_k = cell2mat(H_tilde(k));
                H_bar_m = H_bar_k(:,:,m);
                H_tilde_m = H_tilde_k(:,:,m);
                H_m = H_bar_m + H_tilde_m;
                [~,~,V] = svd(H_m);
                wb= exp( 1j*angle(V(:,1)) ) / sqrt(N_Bm);
                % compute combiner
                for n = 1:N
                    H_bar_m_n = H_bar_k((n-1)*N_Sn+(1:N_Sn),:,m);
                    H_tilde_m_n = H_tilde_k((n-1)*N_Sn+(1:N_Sn),:,m);
                    H_m_n = H_bar_m_n + H_tilde_m_n;
                    wB = wb;
                    ws = exp( -1j * angle(H_m_n * wB) ) / sqrt(N_Sn);
                    SNRm = P*abs(ws.'*H_m_n*wb)^2/sigma2;
                    % compute sum-rate
                    data_sumrate(m) = data_sumrate(m) + df*log2(1+SNRm);
                end 
            end
        end
    end
    % save the BS label and the sum-rate value 
    [cp.m_sumrate, cp.m_tilde] = max(data_sumrate);
    
   
    % generate precoders, index: [AE, BS]
    cp.wb = ones(N_Bm,1)/sqrt(N_Bm);
    k = floor(K/2);
    m = cp.m_tilde;
    if sum(visi_mat(m,:)) > 1e-10
        H_bar_k = cell2mat(H_bar(k));
        H_tilde_k = cell2mat(H_tilde(k));
        H_bar_m = H_bar_k(:,:,m);
        H_tilde_m = H_tilde_k(:,:,m);
        H_m = H_bar_m + H_tilde_m;
        [~,~,V] = svd(H_m);
        cp.wb = exp( 1j*angle(V(:,1)) ) / sqrt(N_Bm);
    end
    % generate combiners, index: [AE, SA]
    cp.ws = ones(N_Sn,N)/sqrt(N_Sn);
    m = cp.m_tilde;
    for n = 1:N  
        if sum(visi_mat(m,n)) > 1e-10
            H_bar_k = cell2mat(H_bar(k));
            H_tilde_k = cell2mat(H_tilde(k));
            H_bar_m_n = H_bar_k((n-1)*N_Sn+(1:N_Sn),:,m);
            H_tilde_m_n = H_tilde_k((n-1)*N_Sn+(1:N_Sn),:,m);
            H_m_n = H_bar_m_n + H_tilde_m_n;
            wB = cp.wb;
            cp.ws(:,n) = exp( -1j * angle(H_m_n * wB) ) / sqrt(N_Sn);
        end
    end
else
    error('Not vaild pcp.func value');
end


end

