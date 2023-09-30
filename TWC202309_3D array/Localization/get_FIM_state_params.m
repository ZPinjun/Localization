function FIM_r = get_FIM_state_params(FIM_eta,sp)

%% get parameters
visi_mat = sp.visi_mat;
Meas = sp.Meas;
p.UEPos = sp.Pu;
p.UEOriR = sp.Ru;

%% calculate transformation matrix
D = sum(sum(visi_mat));
Meas(:,logical(reshape(1-visi_mat.',1,[]))) = [];
T = zeros(5*D,13);
for i = 1:D
    m = Meas(1,i);
    n = Meas(2,i);
    p.D_local = sp.Ps_local(:,n); 
    p.R_local = sp.Rs_local(:,:,n);
    p.BSPos = sp.Pb(:,m);
    p.BSOriR = sp.Rb(:,:,m);
    T_mn = get_trans_mat_mn(p);
    T(([0,D,2*D,3*D,4*D]+i),:) = T_mn;
end

%% calculate FIM of state parameters
FIM_r = T.' * FIM_eta * T;

end

