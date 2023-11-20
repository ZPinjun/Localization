function sp = default_setup()


%% system setup
% --- BS
sp.Pb = [0, 5.01, 3].';
sp.Ob_Euler = [-90 0 0].';
sp.BS_dim = [10,10];
sp.RFC_dim = [5, 5];
sp.BS_spacing = 0.5;   % BS array spacing, in [lambdac]
% --- UE
sp.Pu = [3 2 1].';
% --- RIS
sp.Pr = [-5.01 0 3].';
sp.Or_Euler = [0 0 0].';
sp.RIS_dim = [15, 15].';
sp.RIS_spacing = 0.1;   % RIS array spacing, in [lambdac]
sp.rho = 1e-7;          % clock offset, in [sec]


%% channel setup
sp.fc = 28e9;               % center frequency
sp.BW = 100e6;              % Bandwidth
sp.K = 32;                  % # of subcarriers
sp.G = 9;                   % # of transmissions
sp.P = db2pow(10);          % average transmit power in [mW]
sp.P_RIS = db2pow(7);       % average power supply of active RIS in [mW]
sp.noiseType = 'true';      % noise model: 'true' for correct colored noise, 'aprx' for approximated white Gaussian noise 


%% environment parameters
sp.T = 298.15;               % System temperature (Kelvin), 298.15 Kelvin = 25 celsius
sp.NoiseFigure = 10;         % noise figure in [dB]


%% constants
sp.c = 2.9979e8;             % Speed of light in vacuum
sp.Kb = 1.3806e-23;          % Boltzmann constant

end
