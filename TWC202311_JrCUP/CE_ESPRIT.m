function [parameters_LOS,parameters_RIS] = CE_ESPRIT(y, s, sp)

%% get parameters
N1 = sp.RFC_dim(1); 
N2 = sp.RFC_dim(2);
Nb1 = sp.BS_dim(1);
Nb2 = sp.BS_dim(2);
df = sp.df;
c = sp.c;
Tcell = {sp.T1, sp.T2};
Gamcell = {sp.Gamma1, sp.Gamma2};
BS_spacing = sp.BS_spacing;
RIS_spacing = sp.RIS_spacing;


%% obatain channel matrices
K = size(y,2);
G = size(y,3);
H = zeros(N1*N2,K,G);
for k = 1:K
    for g = 1:G
        H(:,k,g) = y(:,k,g)*conj(s(k,g))/sp.P;
    end
end


%% estimate theta_L_az, theta_L_el, theta_R_az, theta_R_el, tau_L, tau_R
H1 = reshape(H,[N1,N2,K,G]);
H1 = sum(H1, 4);
H1 = permute(H1, [2 1 3]);
Uhat = cpd(H1, 2);
results = zeros(2,3);
for i = 1:3
    if i <= 2
        U1 = Uhat{i};
        T = Tcell{i};
        T1 = T(1:end-1,:);
        T2 = T(2:end,:);
        hatF = pinv(T2)*T1;
        Bh = Tcell{i}';
        tt = [Bh(:,end) hatF'*Bh(:,1)];
        [Qt,~] = qr(tt,0);
        hatQ = eye(size(Tcell{i},2)) - Qt*Qt';
        hatGamma = pinv(hatQ*U1)*hatQ*hatF'*U1;
        e = eig(hatGamma);
        results(:,i) = angle(e);
    else
        U1 = Uhat{i};
        Ur1 = U1(1:end-1,:);
        Ur2 = U1(2:end,:);
        hatPhi = pinv(Ur1)*Ur2;
        e = eig(hatPhi);
        results(:,i) = angle(e);
    end
end
ele = asind(results(:,2)/(BS_spacing*2*pi));
azi = asind(results(:,1)./(BS_spacing*2*pi)./cosd(ele));
delay = -(results(:,3))/(2*pi*df/c);



%% distinguish the channels
T1 = Tcell{1};
T2 = Tcell{2};
W = kron(T1,T2);
a1 = kron( get_a(results(1,1),Nb1), get_a(results(1,2),Nb2) );
a2 = kron( get_a(results(2,1),Nb1), get_a(results(2,2),Nb2) );
aK1 = get_a(results(1,3), K);
aK2 = get_a(results(2,3), K);
AA = [kron(aK1, W'*a1), kron(aK2, W'*a2)];
bb = reshape(H,[N1*N2*K,G]);
alpha1 = zeros(G,1);
alpha2 = zeros(G,1);
for g = 1:G
    alpha_hat = (AA'*AA)^(-1)*AA'*bb(:,g);
    alpha1(g) = alpha_hat(1);
    alpha2(g) = alpha_hat(2);
end
if var(alpha1) < var(alpha2)
    betaR = alpha2;
    indexL = 1;
    indexR = 2;
else 
    betaR = alpha1;
    indexL = 2;
    indexR = 1;
end



%% estimate vartheta_2 & vartheta_3
H1 = reshape(betaR,[sqrt(G) sqrt(G)]).';
Uhat = cpd(H1, 1);
results = zeros(1,2);
for i = 1:2
    U1 = Uhat{i};
    T = Gamcell{i};
    T1 = T(1:end-1,:);
    T2 = T(2:end,:);
    hatF = pinv(T2)*T1;
    Bh = Gamcell{i}';
    tt = [Bh(:,end) hatF'*Bh(:,1)];
    [Qt,~] = qr(tt,0);
    hatQ = eye(size(Gamcell{i},2)) - Qt*Qt';
    hatGamma = pinv(hatQ*U1)*hatQ*hatF'*U1;
    e = eig(hatGamma);
    results(i) = angle(e);
end
vartheta = results/(RIS_spacing*2*pi);

parameters_LOS = [azi(indexL),ele(indexL),delay(indexL)].';
parameters_RIS = [azi(indexR),ele(indexR),delay(indexR),vartheta].';


%% ================================ functions
    function a = get_a(omega,n)
        a = exp(1j*omega*(0:(n-1)).');
    end

end

