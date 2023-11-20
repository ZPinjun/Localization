function cost = CE_cost_fun(X, y, sp)

%% get eta
sp.phiL_az = X(1);
sp.phiL_el = X(2);
sp.phiR_az = X(4);
sp.phiR_el = X(5);
sp.tauL = X(3);
sp.tauR = X(6);
sp.vartheta2 = X(7);
sp.vartheta3 = X(8);


%% get parameters
c = sp.c;
K = sp.K;
G = sp.G;
df = sp.df;
tauR = sp.tauR;
tauL = sp.tauL;
W = sp.W;
lambdac = sp.lambda;
pAE_BS = sp.pAE_BS;
pAE_RIS = sp.pAE_RIS;
phiUB_az = sp.phiL_az;
phiUB_el = sp.phiL_el;
phiRB_az = sp.phiR_az;
phiRB_el = sp.phiR_el;
Gamma = sp.Gamma;
vartheta2 = sp.vartheta2;
vartheta3 = sp.vartheta3;
s = sqrt(sp.P);


%% calculate cost
zL = zeros(size(W,2),K,G);
zR = zeros(size(W,2),K,G);
a_UB_AOA = get_array_response(lambdac, pAE_BS, [phiUB_az, phiUB_el]);
a_RB_AOA = get_array_response(lambdac, pAE_BS, [phiRB_az, phiRB_el]);
tR_tilde = [0, vartheta2, vartheta3].';
PS_RIS_coe = exp((1j*2*pi/lambdac)*tR_tilde.'*pAE_RIS);
for k=1:K
    uL = exp(-1j*2*pi*(k-1)*df*tauL/c)*a_UB_AOA;
    uR_coe = exp(-1j*2*pi*(k-1)*df*tauR/c)*a_RB_AOA;
    zLk = W'*uL;
    for g = 1:G
        % LOS zL
        zL(:,k,g) = zLk;
        % NLOS(RIS) zR
        PS_RIS = PS_RIS_coe*conj(Gamma(:,g));
        uR = uR_coe*PS_RIS;
        zR(:,k,g) = W'*uR;
    end
end

zLv = s*reshape(zL,[],1);
zRv = s*reshape(zR,[],1);
yv = reshape(y,[],1);
% get hat_gL & hat_gR
nzL = norm(zLv,2)^2;
nzR = norm(zRv,2)^2;
nLR = abs(zLv'*zRv)^2;
denominator = nzL*nzR - nLR;
hat_gL = ( zLv'*yv*nzR - (zRv'*yv)*(zLv'*zRv) ) / denominator;
hat_gR = ( zRv'*yv*nzL - (zLv'*yv)*(zRv'*zLv) ) / denominator;

cost = norm( (yv - hat_gL*zLv - hat_gR*zRv), 2)^2;




%% ================================ functions
    function a = get_array_response(lambdac, pE, phi)
        t = get_dir_from_angle(phi);
        a = exp((1j*2*pi/lambdac)*t.'*pE).';
    end



end



