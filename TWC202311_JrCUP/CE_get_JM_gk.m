function D_eta_ch = CE_get_JM_gk(k,g,sp)

%% get parameters
c = sp.c;
df = sp.df;
s = sp.s;
W = sp.W;
Gamma = sp.Gamma;
lambdac = sp.lambda;
pAE_BS = sp.pAE_BS;
pAE_RIS = sp.pAE_RIS;
gL = sp.gL;
gR = sp.gR;
tauR = sp.tauR;
tauL = sp.tauL;
phiUB_az = sp.phiL_az;
phiUB_el = sp.phiL_el;
phiRB_az = sp.phiR_az;
phiRB_el = sp.phiR_el;
vartheta2 = sp.vartheta2;
vartheta3 = sp.vartheta3;


%% compute the Jacobian matrix
x = s(k,g);
[daL_az, daL_el, aL] = get_da(lambdac, pAE_BS, [phiUB_az,phiUB_el]);
[daR_az, daR_el, aR] = get_da(lambdac, pAE_BS, [phiRB_az,phiRB_el]);
[dbeta_dvartheta2,dbeta_dvartheta3,betaR_g,PS_RIS] = get_dvartheta(lambdac,pAE_RIS,vartheta2,vartheta3,Gamma(:,g),gR);

d_thetaL_az = x*W'*gL*exp(-1j*2*pi*(k-1)*df*tauL/c)*daL_az;
d_thetaL_el = x*W'*gL*exp(-1j*2*pi*(k-1)*df*tauL/c)*daL_el;

d_thetaR_az = x*W'*betaR_g*exp(-1j*2*pi*(k-1)*df*tauR/c)*daR_az;
d_thetaR_el = x*W'*betaR_g*exp(-1j*2*pi*(k-1)*df*tauR/c)*daR_el;

d_tauL = x*W'*gL*exp(-1j*2*pi*(k-1)*df*tauL/c)*(-1j*2*pi*(k-1)*df/c)*aL;
d_tauR = x*W'*betaR_g*exp(-1j*2*pi*(k-1)*df*tauR/c)*(-1j*2*pi*(k-1)*df/c)*aR;

d_vartheta2 = x*W'*dbeta_dvartheta2*exp(-1j*2*pi*(k-1)*df*tauR/c)*aR;
d_vartheta3 = x*W'*dbeta_dvartheta3*exp(-1j*2*pi*(k-1)*df*tauR/c)*aR;

d_RgL = x*W'*exp(-1j*2*pi*(k-1)*df*tauL/c)*aL;
d_IgL = 1j*x*W'*exp(-1j*2*pi*(k-1)*df*tauL/c)*aL;
d_RgR = x*W'*exp(-1j*2*pi*(k-1)*df*tauR/c)*aR*PS_RIS;
d_IgR = 1j*x*W'*exp(-1j*2*pi*(k-1)*df*tauR/c)*aR*PS_RIS;

D_eta_ch = [d_thetaL_az, d_thetaL_el, d_thetaR_az, d_thetaR_el, d_tauL, d_tauR, d_vartheta2, d_vartheta3, d_RgL, d_IgL, d_RgR, d_IgR];


%% ================================ functions
    function a = get_array_response(lambdac, pE, phi)
        t = get_dir_from_angle(phi);
        a = exp((1j*2*pi/lambdac)*t.'*pE).';
    end

    function [da_az, da_el, a] = get_da(lambdac, pE, phi)
        a = get_array_response(lambdac, pE, phi);
        t1 = [0, cosd(phi(1))*cosd(phi(2)), 0];
        t2 = [0, -sind(phi(1))*sind(phi(2)), cosd(phi(2))];
        da_az = a .* ((1j*2*pi/lambdac)*t1*pE).';
        da_el = a .* ((1j*2*pi/lambdac)*t2*pE).';
    end

    function [dbeta_dvartheta2,dbeta_dvartheta3,betaR_g,PS_RIS] = get_dvartheta(lambdac,pE,v2,v3,gamma,gR)
        tR_tilde = [0, v2, v3].';
        PS_RIS = exp((1j*2*pi/lambdac)*tR_tilde.'*pE)*conj(gamma);
        betaR_g = gR*PS_RIS;
        dbeta_dvartheta2 = ( ((1j*2*pi/lambdac)*pE(2,:)).*(exp((1j*2*pi/lambdac)*tR_tilde.'*pE)) ) * conj(gamma) * gR;
        dbeta_dvartheta3 = ( ((1j*2*pi/lambdac)*pE(3,:)).*(exp((1j*2*pi/lambdac)*tR_tilde.'*pE)) ) * conj(gamma) * gR;
    end

end

