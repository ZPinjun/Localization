clc; clear; close all

sp = default_setup();
sp = update_setup(sp);
[y0, s, H, CovMat, H_R2, PowNoise] = gen_channel(sp);
noise = gen_noise(sp, 1, 1, H_R2, PowNoise);
y = y0 + noise(:,:,:,1);

% true channel parameters & localization parameters
tp_ch = [sp.phiL_az,sp.phiL_el,sp.tauL,sp.phiR_az,sp.phiR_el,sp.tauR,sp.vartheta2,sp.vartheta3];
tp_loc = [sp.Pu; sp.Pr; sp.Or_Euler(1); sp.rho*sp.c].';

disp('----------------------- CRLBs -----------------------');
% compute CRLB
[CRLB_eta, CRLB_r] = get_CRLB(sp,1,CovMat);
disp('CRLB of channel parameters:');
disp(['EB(thetaL): ', num2str(CRLB_eta(1)),'; EB(thetaR): ', num2str(CRLB_eta(2)),'; EB(tauL): ', num2str(CRLB_eta(3)), ...
    '; EB(tauR): ', num2str(CRLB_eta(4)),'; EB(vartheta): ', num2str(CRLB_eta(5)), newline]);

disp('CRLB of channel parameters:');
disp(['EB(pU): ', num2str(CRLB_r(1)),'; EB(pR): ', num2str(CRLB_r(2)),'; EB(o3): ', num2str(CRLB_r(3)), ...
    '; EB(Delta): ', num2str(CRLB_r(4)), newline]);


disp('----------------------- Channel estimation -----------------------');
% tensor-ESPRIT
disp('Coarse ==> tensor-ESPRIT:');
[pL, pR] = CE_ESPRIT(y, s, sp);
disp(['thetaL_az: ', num2str(pL(1)),'; thetaL_el: ', num2str(pL(2)),'; tauL: ', num2str(pL(3)), ...
    '; thetaR_az: ', num2str(pR(1)),'; thetaR_el: ', num2str(pR(2)),'; tauR: ', num2str(pR(3)), ...
    '; vartheta_az: ', num2str(pR(4)),'; vartheta_el: ', num2str(pR(5)), newline]);

% SOMP
disp('Coarse ==> SOMP:');
es_SOMP = CE_SOMP(y, sp, 128);   % D = 128^2
disp(['thetaL_az: ', num2str(es_SOMP(1)),'; thetaL_el: ', num2str(es_SOMP(2)),'; tauL: ', num2str(es_SOMP(3)), ...
    '; thetaR_az: ', num2str(es_SOMP(4)),'; thetaR_el: ', num2str(es_SOMP(5)),'; tauR: ', num2str(es_SOMP(6)), ...
    '; vartheta_az: ', num2str(es_SOMP(7)),'; vartheta_el: ', num2str(es_SOMP(8)), newline]);

% tensor-ESPRIT + LS refinement
disp('Refined ==> tensor-ESPRIT + LS:');
Xinit = [pL;pR];
X_ESPRIT = CE_refinement(Xinit, y, sp);
disp(['thetaL_az: ', num2str(X_ESPRIT(1)),'; thetaL_el: ', num2str(X_ESPRIT(2)),'; tauL: ', num2str(X_ESPRIT(3)), ...
    '; thetaR_az: ', num2str(X_ESPRIT(4)),'; thetaR_el: ', num2str(X_ESPRIT(5)),'; tauR: ', num2str(X_ESPRIT(6)), ...
    '; vartheta_az: ', num2str(X_ESPRIT(7)),'; vartheta_el: ', num2str(X_ESPRIT(8)), newline]);

% SOMP + LS refinement
disp('Refined ==> SOMP + LS:');
Xinit = es_SOMP.';
X_SOMP = CE_refinement(Xinit, y, sp);
disp(['thetaL_az: ', num2str(X_SOMP(1)),'; thetaL_el: ', num2str(X_SOMP(2)),'; tauL: ', num2str(X_SOMP(3)), ...
    '; thetaR_az: ', num2str(X_SOMP(4)),'; thetaR_el: ', num2str(X_SOMP(5)),'; tauR: ', num2str(X_SOMP(6)), ...
    '; vartheta_az: ', num2str(X_SOMP(7)),'; vartheta_el: ', num2str(X_SOMP(8)), newline]);

% true channel parameters
disp('True channel parameters:');
disp(['thetaL_az: ', num2str(tp_ch(1)),'; thetaL_el: ', num2str(tp_ch(2)),'; tauL: ', num2str(tp_ch(3)), ...
    '; thetaR_az: ', num2str(tp_ch(4)),'; thetaR_el: ', num2str(tp_ch(5)),'; tauR: ', num2str(tp_ch(6)), ...
    '; vartheta_az: ', num2str(tp_ch(7)),'; vartheta_el: ', num2str(tp_ch(8)), newline]);

% Note: the units of angles and delays are degree and meter, respectively

disp('----------------------- Conversion to localization domain -----------------------');
% 2D search
disp('Conversion to localization domain:');
% --------------------------------------------
% these search range should be properly set!!!
olim = [-30, 30];
dlim = [28,32];
% --------------------------------------------
r = Loc_coarse(X_ESPRIT.', sp, olim, dlim);   % input the results of tensor-ESPRIT + LS refinement
r = r.';
disp(['pu: [', num2str(r(1:3)),']; pr: [', num2str(r(4:6)),']; o3: ', num2str(r(7)), '; Delta: ', num2str(r(8)), newline]);

disp('True localization parameters:');
disp(['pu: [', num2str(tp_loc(1:3)),']; pr: [', num2str(tp_loc(4:6)),']; o3: ', num2str(tp_loc(7)), '; Delta: ', num2str(tp_loc(8))]);

% Note: the position {pU, pR} in units of meter, the Euler angle o3 in units of degree, and the clock bias Delta in units of meter
