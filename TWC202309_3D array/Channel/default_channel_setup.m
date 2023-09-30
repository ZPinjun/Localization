% function to generate the default channel parameters

function cp = default_channel_setup()

% constants
cp.c = 2.9979e8;             % Speed of light in vacuum
cp.h = 6.6262e-34;           % Planck constant
cp.Kb = 1.3806e-23;          % Boltzmann constant
cp.R = 8.3144;               % Gas constant
cp.Na = 6.0221e23;           % Avogadro constant
cp.Tstp = 273.15;            % Temperature at standard pressure

cp.Kr = 4;   % K-factor of Ricean fading
cp.fc = 140e9;   % center frequancy
cp.BW = 1000e6;   % Bandwidth
cp.K = 128;   % # of subcarriers
cp.G = 10;   % # of transmissions
cp.nu = 2;   % path loss component
cp.P = 10;   % average transmit power in [mW]

% Molecules of transmission medium
cp.molecules = {'N2','O2','H2O','CO2','CH4'};
% Molecules ratio
cp.moleculesRatio = [76.6 21.0 1.6 0.03 0.77]/100;

cp.T0 = 296;                 % Reference temperature (Kelvin)
cp.T = 298.15;               % System temperature (Kelvin), 298.15 Kelvin = 25 celsius
cp.p0 = 1;                   % Standard pressure (atm)
cp.p = cp.p0;                % System pressure (atm)


end

