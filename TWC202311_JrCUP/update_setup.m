function sp = update_setup(sp)

rng(1)

%% get parameters
% --- BS parameters
Pb = sp.Pb;                         
Ob_Euler = sp.Ob_Euler;             
BS_dim = sp.BS_dim;
BS_spacing = sp.BS_spacing;
RFC_dim = sp.RFC_dim;
% --- UE parameters
Pu = sp.Pu;  
% --- RIS parameters
Pr = sp.Pr;                         
P_RIS = sp.P_RIS;
RIS_dim = sp.RIS_dim;
RIS_spacing = sp.RIS_spacing;
Or_Euler = sp.Or_Euler;  
% --- channel parameters
rho = sp.rho;                       
c = sp.c;
P = sp.P;
fc = sp.fc;
G = sp.G;
K = sp.K;


%% generate BS & RIS orientations
Rb = eul2rotm(reshape(deg2rad(Ob_Euler), [1,3]), 'ZYX');
Rr = eul2rotm(reshape(deg2rad(Or_Euler), [1,3]), 'ZYX');

%% generate channel parameters:
% vector structure: [tauL, tauR, phiL_az, phiL_el, phiR_az, phiR_el, vartheta2, vartheta3, gL, gR]
dbu = norm(Pb-Pu,2);
dbr = norm(Pb-Pr,2);
dur = norm(Pu-Pr,2);
lambda = c/fc;
% LOS delay, in [m]
tauL = dbu + rho*c;  
% NLOS delay, in [m]
tauR = dbr + dur + rho*c;  
% AOA from UE to BS
tBU = Rb.'*( (Pu - Pb)/dbu);
[phiL_az, phiL_el] = get_angle_from_dir(tBU); 
% AOA from RIS to BS
tBR = Rb.'*( (Pr - Pb)/dbr);
[phiR_az, phiR_el] = get_angle_from_dir(tBR);   
% intermediate unknowns
tRU = Rr.'*( (Pu - Pr)/dur);
tRB = Rr.'*( (Pb - Pr)/dbr);
vartheta2 = tRU(2) + tRB(2);
vartheta3 = tRU(3) + tRB(3);
% LOS channel gain
gL = lambda*exp(1j*rand(1)*2*pi) / (4*pi*dbu);   
% NLOS channel gain
gR1 = lambda/(4*pi*dur);  
gR2 = lambda/(4*pi*dbr);  
gR = gR1*gR2*exp(1j*rand(1)*2*pi);


%% generate the position of BS elements
% arrays are deployed on the Y-O-Z plane, x-axis is the normal direction
yrange = (  (1:BS_dim(1))-1  ) * (BS_spacing*lambda);
zrange = (  (1:BS_dim(2))-1  ) * (BS_spacing*lambda);
pAE_BS = [zeros(1,BS_dim(1)*BS_dim(2));
        kron(yrange, ones(1,BS_dim(1)));
        kron(ones(1,BS_dim(1)), zrange)];


%% generate the position of RIS elements
% arrays are deployed on the Y-O-Z plane, x-axis is the normal direction
yrange = (  (1:RIS_dim(1)) - 1  ) * (RIS_spacing*lambda);
zrange = (  (1:RIS_dim(2)) - 1  ) * (RIS_spacing*lambda);
pAE_RIS = [zeros(1,RIS_dim(1)*RIS_dim(2));
        kron(yrange, ones(1,RIS_dim(1)));
        kron(ones(1,RIS_dim(1)), zrange)];

    
%% generate noise
% thermal noise at receiver
sp.N0 = sp.Kb*sp.T*1000;                                % thermal noise PSD in [mW/Hz]
sp.Pn = sp.N0*sp.BW;                                    % thermal noise in [mW]
sp.sigma_in = sqrt(sp.Pn);                              % input noise sigma
sp.sigma = sqrt(10^(sp.NoiseFigure/10))*sp.sigma_in;    % thermal noise at the receiver
% thermal noise at active RIS
sp.sigma_r = sp.sigma;


%% generate RIS profile & combiner 
% RIS
RIS_amp =sqrt( (P_RIS) / (RIS_dim(1)*RIS_dim(2)*( P*gR1^2 + sp.sigma_r^2 ) ) + 1 );
gamma1 = rand(1,sqrt(G));
gamma2 = rand(1,sqrt(G));
Gamma1 = exp(1j*(1:RIS_dim(1)).'*gamma1);
Gamma2 = exp(1j*(1:RIS_dim(2)).'*gamma2);
Gamma = RIS_amp*kron(Gamma1, Gamma2);
% combiner, satisfying shift invariance property
w1 = randn(1,RFC_dim(1));
w2 = randn(1,RFC_dim(2));
T1 = exp(1j*(1:BS_dim(1)).'*w1);
T2 = exp(1j*(1:BS_dim(2)).'*w2);
W = kron(T1, T2);
% generate a random RIS profile & combiner (for OMP-based channel estimation)
Gamma_random = RIS_amp*exp(1j*randn(size(Gamma)));
W_random = exp(1j*randn(size(W)));

%% generate transmitted signals
sp.s = sqrt(P)*ones(K,G);


%% get frequency fk, k = 1,...,K
df = sp.BW/sp.K;                            % Subcarrier bandwidth (Hz)
fstart_sub = sp.fc - sp.BW/2 + df/2;
fstop_sub = fstart_sub + sp.BW - df/2;
sp.fk = fstart_sub:df:fstop_sub;            % Center frequency of each subcarrier (Hz)
sp.df = df;



%% update setup
% --- BS 
sp.Rb = Rb;
sp.pAE_BS = pAE_BS;
% --- RIS
sp.Rr = Rr;
sp.pAE_RIS = pAE_RIS;
sp.lambda = lambda;
% --- channel parameters
sp.tauL = tauL;
sp.tauR = tauR;
sp.phiL_az = phiL_az;
sp.phiL_el = phiL_el;
sp.phiR_az = phiR_az;
sp.phiR_el = phiR_el;
sp.vartheta2 = vartheta2;
sp.vartheta3 = vartheta3;
sp.gL = gL;
sp.gR = gR;
sp.gR1 = gR1;
sp.gR2 = gR2;
% RIS profile & combiner
sp.T1 = T1;
sp.T2 = T2;
sp.W = W;
sp.W_random = W_random;
sp.Gamma1 = Gamma1;
sp.Gamma2 = Gamma2;
sp.Gamma = Gamma;
sp.RIS_amp = RIS_amp;
sp.Gamma_random = Gamma_random;

end



