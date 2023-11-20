function noise = gen_noise(sp, repeatNum, scale, H_R2, PowNoise)

rng(0)

%% get parameters
c = sp.c;
K = sp.K;
G = sp.G;
N1 = sp.RFC_dim(1);
N2 = sp.RFC_dim(2);

if strcmp(sp.noiseType, 'true')
    % generate correct colored noise
    tauR = sp.tauR;
    df = sp.df;
    NB = size(H_R2, 1);
    NR = size(H_R2, 2);
    sigma0 = sp.sigma;
    sigmar = sp.sigma_r;
    W = sp.W;
    Gamma = sp.Gamma;
    noise = zeros(N1*N2,repeatNum,K,G);
    for k=1:K
        Hk_R2 = H_R2*exp(-1j*2*pi*(k-1)*df*tauR/c);
        for g = 1:G
            n0 = (sqrt(0.5)*sigma0*randn(NB, repeatNum) + 1j*sqrt(0.5)*sigma0*randn(NB, repeatNum)) * scale;
            nr = (sqrt(0.5)*sigmar*randn(NR, repeatNum) + 1j*sqrt(0.5)*sigmar*randn(NR, repeatNum)) * scale;
            noise(:,:,k,g) = W'*( Hk_R2*diag(Gamma(:,g))*nr + n0 );
        end
    end
    noise = permute(noise, [1, 3, 4, 2]);
elseif strcmp(sp.noiseType, 'aprx')
    % generate approximated AWGN
    p = PowNoise/(K*G*N1*N2);
    noise = (sqrt(2)/2)*scale*sqrt(p)*(  randn(N1*N2,K,G,repeatNum) + 1j*randn(N1*N2,K,G,repeatNum)  );
else
    error('Unrecognized noise type!')
end

end


