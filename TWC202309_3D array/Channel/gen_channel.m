% function to generate the LOS component of THz channel

function cp = gen_channel(sp,cp)


%% get parameters
M = size(sp.Pb,2);   % # of BSs
N = size(sp.Ps_local,2);   % # of SAs
Pb = sp.Pb;
Rb = sp.Rb;
Ps = sp.Ps;
Rs = sp.Rs;
N_Bm = sp.NB_dim(1)*sp.NB_dim(2);   % # of antenna elements in each BS
N_Sn = sp.NS_dim(1)*sp.NS_dim(2);   % # of antenna elements in each SA
Meas = sp.Meas;
Gb = sp.Gb;
Gs = sp.Gs;
visi_mat = sp.visi_mat;
pAE_B = sp.pAE_B;
pAE_S = sp.pAE_S;

K = cp.K;   % # of subcarriers
Kr = cp.Kr;   % K-factor of Ricean fading
nu = cp.nu;   % path loss component
c = cp.c;
fk = cp.fk;



%% get molecular absorption coefficient for differenct fk
p.fk = fk;
p.K = K;
p.c = c;
p.molecules = cp.molecules;
p.moleculesRatio = cp.moleculesRatio;
p.h = cp.h;           % Planck constant
p.Kb = cp.Kb;         % Boltzmann constant
p.R = cp.R;           % Gas constant
p.Na = cp.Na;         % Avogadro constant
p.Tstp = cp.Tstp;     % Temperature at standard pressure
p.T0 = cp.T0;         % Reference temperature (Kelvin)
p.T = cp.T;           % System temperature (Kelvin) 
p.p0 = cp.p0;         % Standard pressure (atm)
p.p = cp.p;           % System pressure (atm)
K_abs = Abs_Coef_Hitran(p);


%% get phases of the beamsteering vector in [m]
PsiB = zeros(N_Bm,M*N);
PsiS = zeros(N_Sn,M*N);
for m = 1:M
    for n = 1:N
        if visi_mat(m,n) == 1
            % AOD
            tAOD_global = Ps(:,n) - Pb(:,m);   % in the global coordinate system
            tAOD_local = Rb(:,:,m).'*tAOD_global;   % in the BS's local coordinate system
            %tAOD_local = normalize(tAOD_local, 'norm', 2);
            tAOD_local = tAOD_local/norm(tAOD_local,2);
            PsiB(:,(m-1)*N+n) = pAE_B.' * tAOD_local;
            % AOA
            tAOA_global = Pb(:,m) - Ps(:,n);   % in the global coordinate system
            tAOA_local = Rs(:,:,n).'*tAOA_global;   % in the BS's local coordinate system
            %tAOA_local = normalize(tAOA_local, 'norm', 2);
            tAOA_local = tAOA_local/norm(tAOA_local,2);
            PsiS(:,(m-1)*N+n) = pAE_S.' * tAOA_local;
        end
    end
end


%% calculate the LOS component of the channel matrix
H_bar = cell(1,K);
Gmn = zeros(M,N,K);   % save G_mn values
% loop over subcarries
for k = 1:K
    f = fk(k);
    H_bar_k = zeros(N*N_Sn, N_Bm, M);
    for m = 1:M
        for n = 1:N
            G_mn = 0; 
            H_bar_k_mn = zeros(N_Sn, N_Bm);            
            if visi_mat(m,n) == 1
                d_mn = norm(Pb(:,m)-Ps(:,n),2);   % propagation distance
                G_mn = (c/(4*pi*f*d_mn))^(nu/2)*exp(-0.5*K_abs(k)*d_mn)*Gb(m,n)*Gs(m,n);
                tau_mn = Meas(end,(m-1)*N+n)/c;   % delay 
                aB_m = exp(1j*2*pi*f/c*PsiB(:,(m-1)*N+n));
                aS_n = exp(1j*2*pi*f/c*PsiS(:,(m-1)*N+n)); 
                H_bar_k_mn = G_mn * exp(-1j*2*pi*f*tau_mn) * aS_n * aB_m.';
            end
            Gmn(m,n,k) = G_mn;
            H_bar_k((n-1)*N_Sn+(1:N_Sn),:,m) = H_bar_k_mn;
        end
    end
    H_bar(k) = mat2cell(H_bar_k * sqrt(Kr/(Kr+1)), N_Sn*N, N_Bm, M);
end
 

%% calculate the NLOS component of the channel matrix
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



%% channel outputs
cp.H_bar = H_bar;
cp.Gmn = Gmn;
cp.H_tilde = H_tilde;



end

