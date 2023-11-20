function [y, s, H, CovMat, H_R2, PowNoise] = gen_channel(sp)

%% get parameters
K = sp.K;
G = sp.G;
N1 = sp.RFC_dim(1);
N2 = sp.RFC_dim(2);
df = sp.df;
s = sp.s;
W = sp.W;
Gamma = sp.Gamma;
lambdac = sp.lambda;
pAE_BS = sp.pAE_BS;
pAE_RIS = sp.pAE_RIS;
gL = sp.gL;
gR = sp.gR;
gR2 = sp.gR2;
tauR = sp.tauR;
tauL = sp.tauL;
phiUB_az = sp.phiL_az;
phiUB_el = sp.phiL_el;
phiRB_az = sp.phiR_az;
phiRB_el = sp.phiR_el;
vartheta2 = sp.vartheta2;
vartheta3 = sp.vartheta3;
sigma02 = sp.sigma^2;
sigmar2 = sp.sigma_r^2;
c = sp.c;


%% generate channel matrix, received signals, noise covariance matrix
H = zeros(N1*N2,K,G);
y = zeros(N1*N2,K,G);
CovMat = zeros(2*N1*N2, 2*N1*N2, K, G);
a_UB_AOA = get_array_response(lambdac, pAE_BS, [phiUB_az,phiUB_el]);
a_RB_AOA = get_array_response(lambdac, pAE_BS, [phiRB_az,phiRB_el]);
tR_tilde = [0, vartheta2, vartheta3].';
% ----------------------------------------------------
% avoid repeating calculations in loops
InPhase = exp((1j*2*pi/lambdac)*tR_tilde.'*pAE_RIS);
C0 = get_Cn(W', sigma02);
Hkg_LOS1 = gL*W'*a_UB_AOA;
H_R2 = gR2*a_RB_AOA*InPhase;
RISpath_coe1 = gR*W'*a_RB_AOA;
% ----------------------------------------------------
PowNoise = 0;
for k=1:K
    % LOS
    Hkg_LOS = Hkg_LOS1*exp(-1j*2*pi*(k-1)*df*tauL/c);
    % RIS path coefficient
    Hk_R2 = H_R2*exp(-1j*2*pi*(k-1)*df*tauR/c);
    RISpath_coe = RISpath_coe1*exp(-1j*2*pi*(k-1)*df*tauR/c);
    for g = 1:G
        % RIS
        PS_RIS = InPhase*conj(Gamma(:,g));
        Hkg_RIS = RISpath_coe*PS_RIS;
        % complete channel
        Hkg = Hkg_LOS + Hkg_RIS;
        H(:,k,g) = Hkg;
        y(:,k,g) = Hkg*s(k,g);
        % --- compute noise power for Gaussian approximation
        W1 = W'*Hk_R2*diag(Gamma(:,g));
        CovMat(:,:,k,g) = C0 + get_Cn(W1, sigmar2); 
        PowNoise = PowNoise + trace(CovMat(:,:,k,g));
    end
end



%% ================================ functions
    function a = get_array_response(lambdac, pE, phi)
        t = get_dir_from_angle(phi);
        a = exp((1j*2*pi/lambdac)*t.'*pE).';
    end

    function Cn = get_Cn(A, sigma2)
        AR = real(A);
        AI = imag(A);
        Cn = 0.5*sigma2*[AR*AR.'+ AI*AI.', AR*AI.' - AI*AR.'; AI*AR.' - AR*AI.', AR*AR.'+ AI*AI.'];
    end


end

