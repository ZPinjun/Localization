% function to generate the default system parameters

function sp = default_system_setup()


sp.fc = 140e9;   % default center frequency for system in [Hz]
sp.c = 2.9979e8;  % Speed of light in vacuum

% generate BS
sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5].';
sp.Ob_Euler = [45 135 0; 90 0 0].';

% generate UE
sp.Pu = [0 0 0].';
sp.Ou_Euler = [0 0 0].';

% type of subarray
sp.TypeSA = 'cuboidal';

% array size
sp.NB_dim = [10, 10].';
sp.NS_dim = [4, 4].';

sp.rho = 1e-7;   % clock offset, in [sec]

sp.vartheta = 180;   % antenna directivity parameter of cone model (Eq. (8), [1]), in [deg], 180 means a semi-sphere radiation pattern.


% Reference: 
% [1] V. Petrov, M. Komarov, D. Moltchanov, J. M. Jornet and Y. Koucheryavy, "Interference and SINR in Millimeter Wave and Terahertz Communication Systems With Blocking and Directional Antennas," in IEEE Transactions on Wireless Communications, vol. 16, no. 3, pp. 1791-1808, March 2017, doi: 10.1109/TWC.2017.2654351.

end

