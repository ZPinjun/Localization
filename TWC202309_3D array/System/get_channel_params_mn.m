% get channel parameters (ADO, AOA, delay) and visibility of a single path
% MeaPath_mn = [azi_AOD, ele_AOD, azi_AOA, ele_AOA, delay].';

function [MeaPath_mn, visibility] = get_channel_params_mn(PBS, RBS, PUE, RUE, rho, c, vartheta)

% AOD direction vector
tBU = PUE - PBS;   % in the global coordinate system
tBU_BSCS = RBS.'*tBU;   % in the BS's local coordinate system 
%t_AOD = normalize(tBU_BSCS, 'norm', 2);
t_AOD = tBU_BSCS/norm(tBU_BSCS,2);

% AOA direction vector
tUB = PBS - PUE;   % in the global coordinate system
tUB_UECS = RUE.'*tUB;   % in the UE's local coordinate system  
%t_AOA = normalize(tUB_UECS, 'norm', 2);
t_AOA = tUB_UECS/norm(tUB_UECS,2);

% get channel parameters
if acosd(t_AOD(1)) < (vartheta/2 - 1e-10) && acosd(t_AOA(1)) < (vartheta/2 - 1e-10)
    visibility = 1;
    AOD = get_angle_from_dir(tBU_BSCS);
    AOA = get_angle_from_dir(tUB_UECS);
    tau = norm(tUB, 2) + rho*c;  % delay in [m]
    MeaPath_mn = [AOD; AOA; tau];    
else 
    visibility = 0;
    MeaPath_mn = nan(5,1);
end












end

