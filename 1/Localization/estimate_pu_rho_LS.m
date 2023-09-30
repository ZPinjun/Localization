% Meas: angle and delay measurements
% Ru: UE orientation
% np: known parameters
% Ues unit [m] for delay tau and clock offset rho

function pu_rho_est = estimate_pu_rho_LS(Ru, np)

Meas = np.Meas;
Q = size(Meas,2);
Rbs = np.Rbs;
Rsa = np.Rsa;
pbs = np.pbs;
dsa = np.dsa;


A = zeros(Q*3, 4);
B = zeros(Q*3, 4);
C = zeros(Q*3, 1);
D = zeros(Q*3, 1);
for i = 1:Q
    Rbm = Rbs(:,:,Meas(1,i));
    Rn = Rsa(:,:,Meas(2,i));
    a = Rbm*[cos(Meas(2+2,i))*cos(Meas(1+2,i)); cos(Meas(2+2,i))*sin(Meas(1+2,i)); sin(Meas(2+2,i))];
    a1 = -Ru*Rn*[cos(Meas(4+2,i))*cos(Meas(3+2,i)); cos(Meas(4+2,i))*sin(Meas(3+2,i)); sin(Meas(4+2,i))];
    A((1:3)+3*(i-1),:) = [eye(3), a];
    B((1:3)+3*(i-1),:) = [eye(3), a1];

    pbm = pbs(:,Meas(1,i));
    dn = dsa(:,Meas(2,i));
    tau = Meas(5+2, i);   % in [m]
    C((1:3)+3*(i-1),:) = pbm + a*tau - Ru*dn;
    D((1:3)+3*(i-1),:) = pbm + a1*tau - Ru*dn;
end
    G = [A;B];
    b = [C;D];
    pu_rho_est = (G.'*G)^(-1)*G.'*b;
end

