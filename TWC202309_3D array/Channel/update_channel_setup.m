function cp = update_channel_setup(cp)

%% generate noise
cp.N0 = cp.Kb*cp.T*1000;   % thermal noise PSD in [mW/Hz]
cp.operationBW = cp.BW;    % Operation bandwidth for Thermal noise
cp.Pn = cp.N0*cp.operationBW;    % thermal noise in [mW]
cp.Pn_dBm = 10*log10(cp.Pn);    % thermal noise decibel in [dBm]
cp.sigma_in = sqrt(cp.Pn);    % input noise sigma
cp.NoiseFigure = 10;       % noise figure in [dB].
cp.sigma = sqrt(10^(cp.NoiseFigure/10))*cp.sigma_in;

%% get frequency fk, k = 1,...,K
df = cp.BW/cp.K;   % Subcarrier bandwidth (Hz)
fstart_sub = cp.fc - cp.BW/2 + df/2;
fstop_sub = fstart_sub + cp.BW - df/2;
cp.fk = fstart_sub:df:fstop_sub;  % Center frequency of each subcarrier (Hz)


end

