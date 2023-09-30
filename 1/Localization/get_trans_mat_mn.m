% get transformation matrix of path m-n

function T_mn = get_trans_mat_mn(p)


%% get parameters
% UE
pu = p.UEPos;
Ru = p.UEOriR;
% SA
dn = p.D_local;
Rn = p.R_local;
% BS
pBm = p.BSPos;
RBm = p.BSOriR;

u1 = [1 0 0].';
u2 = [0 1 0].';
u3 = [0 0 1].';
tbs = pu + Ru*dn - pBm;
dir_BS = RBm.'*tbs;
dir_SA = Rn.'*Ru.'*tbs;
normtbs = norm(tbs,2);
termdn = [dn(1)*eye(3); dn(2)*eye(3); dn(3)*eye(3)];


%% row 1 (AOD_azimuth)
dAODazi_dpu = (dir_BS(1)*RBm*u2 - dir_BS(2)*RBm*u1) / (dir_BS(1)^2+dir_BS(2)^2);
dAODazi_dRu = (dir_BS(1)*RBm*u2*dn.' - dir_BS(2)*RBm*u1*dn.') / (dir_BS(1)^2+dir_BS(2)^2);
dAODazi_dRu_vec = reshape(dAODazi_dRu, [9,1]);
row1 = [dAODazi_dpu.', dAODazi_dRu_vec.', 0];

%% row 2 (AOD_elevation)
term2 = 1/sqrt(1 - (dir_BS(3)/normtbs)^2);
dAODele_dpu = term2 * (  RBm*u3/normtbs - dir_BS(3)*tbs/normtbs^3  );
dAODele_dRu_vec = term2 * (  reshape(RBm*u3*dn.', [9,1])/normtbs - dir_BS(3)*termdn*tbs/normtbs^3  );
row2 = [dAODele_dpu.', dAODele_dRu_vec.', 0];

%% row3 (AOA_azimuth)
dAOAazi_dpu = (dir_SA(1)*Ru*Rn*u2 - dir_SA(2)*Ru*Rn*u1)/(dir_SA(1)^2+dir_SA(2)^2);
dtermu1 =  (pu - pBm)*u2.'*Rn.' + Ru*(Rn*u2*dn.' + dn*u2.'*Rn.');
dtermu2 =  (pu - pBm)*u1.'*Rn.' + Ru*(Rn*u1*dn.' + dn*u1.'*Rn.');
dAOAazi_dRu = (dir_SA(1)*dtermu1 - dir_SA(2)*dtermu2)/(dir_SA(1)^2+dir_SA(2)^2);
dAOAazi_dRu_vec = reshape(dAOAazi_dRu, [9,1]);
row3 = [dAOAazi_dpu.', dAOAazi_dRu_vec.', 0];

%% row4 (AOA_elevation)
term4 = 1/sqrt(1 - (dir_SA(3)/normtbs)^2);
dAOAele_dpu = -term4 * ( Ru*Rn*u3/normtbs - dir_SA(3)*tbs/normtbs^3 );
dtermu =  (pu - pBm)*u3.'*Rn.' + Ru*(Rn*u3*dn.' + dn*u3.'*Rn.');
dAOAele_dRu_vec = -term4 * (   reshape(dtermu, [9,1])/normtbs - dir_SA(3)*termdn*tbs/normtbs^3   );
row4 = [dAOAele_dpu.', dAOAele_dRu_vec.', 0];

%% row5 (delay)
dtau_dpu = tbs/normtbs;  
dtau_dRu_vec = termdn*tbs/normtbs;
row5 = [dtau_dpu.', dtau_dRu_vec.', 1];

T_mn = [row1; row2; row3; row4; row5];

end

