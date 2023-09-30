function sp = update_system_setup(sp)

%% generate rotation matrix
% BSs
sp.Rb = zeros(3,3,size(sp.Pb,2));
for m = 1:size(sp.Pb,2)
    sp.Rb(:,:,m) = eul2rotm(reshape(deg2rad(sp.Ob_Euler(:,m)),[1,3]), 'ZYX');
end
% UE
sp.Ru = eul2rotm(reshape(deg2rad(sp.Ou_Euler), [1,3]), 'ZYX');


%% generate SAs
[Dsubarray,Rsubarray] = gen_SAs(sp.TypeSA);
sp.Ps_local = Dsubarray;   % in the UE's local coordinate system
sp.Rs_local = Rsubarray;   % in the UE's local coordinate system
M = size(sp.Pb,2);
N = size(Dsubarray,2);
sp.Ps = sp.Pu + sp.Ru*Dsubarray;
sp.Rs = zeros(3,3,N);
for n = 1:N
    sp.Rs(:,:,n) = sp.Ru*Rsubarray(:,:,n);
end


%% generate channel parameters: AOD, AOA, delay
Meas = zeros(7, M*N);   % row 3~7: AOD_azi, AOd_ele, AOA_azi, AOA_ele, delay
Meas(1,:) = kron(1:M,ones(1,N));   % the first row: the label of BS
Meas(2,:) = kron(ones(1,M),1:N);   % the second row: the label of SA
visib_mat = zeros(M,N);   % matrix indicating the visibility of the paths
for m = 1:M
    PBS = sp.Pb(:,m);
    RBS = sp.Rb(:,:,m);
    for n = 1:N
        PSA = sp.Ps(:,n);
        RSA = sp.Rs(:,:,n);
        [Meas(3:end, (m-1)*N+n), visibility] = get_channel_params_mn(PBS, RBS, PSA, RSA, sp.rho, sp.c, sp.vartheta);
        visib_mat(m,n) = visibility;        
    end
end
sp.Meas = Meas;
sp.visi_mat = visib_mat;


%% generate antenna gains, Eq. (8) in [1]
G = 2 / ( 1 - cosd(sp.vartheta/2) );
sp.Gb = sqrt(G) * visib_mat;
sp.Gs = sqrt(G) * visib_mat;



%% generate the position of antenna elements
% arrays are deployed on the Y-O-Z plane, x-axis is the normal direction
yrange_BS = (  (1:sp.NB_dim(1)) - (1+sp.NB_dim(1))/2  ) * (0.5*sp.c/sp.fc);
zrange_BS = (  (1:sp.NB_dim(2)) - (1+sp.NB_dim(2))/2  ) * (0.5*sp.c/sp.fc);
sp.pAE_B = [zeros(1,sp.NB_dim(1)*sp.NB_dim(2));
        kron(yrange_BS, ones(1,sp.NB_dim(1)));
        kron(ones(1,sp.NB_dim(1)), zrange_BS)];
yrange_SA = (  (1:sp.NS_dim(1)) - (1+sp.NS_dim(1))/2  ) * (0.5*sp.c/sp.fc);
zrange_SA = (  (1:sp.NS_dim(2)) - (1+sp.NS_dim(2))/2  ) * (0.5*sp.c/sp.fc);
sp.pAE_S = [zeros(1,sp.NS_dim(1)*sp.NS_dim(2));
        kron(yrange_SA, ones(1,sp.NS_dim(1)));
        kron(ones(1,sp.NS_dim(1)), zrange_SA)];


% Reference: 
% [1] V. Petrov, M. Komarov, D. Moltchanov, J. M. Jornet and Y. Koucheryavy, "Interference and SINR in Millimeter Wave and Terahertz Communication Systems With Blocking and Directional Antennas," in IEEE Transactions on Wireless Communications, vol. 16, no. 3, pp. 1791-1808, March 2017, doi: 10.1109/TWC.2017.2654351.

end

