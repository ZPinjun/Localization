function [r_est, cost] = Loc_coarse(eta, sp, olim, dlim, oNum, dNum)

% get parameters
pb = sp.Pb;
Rb = sp.Rb;
azL = eta(1);
elL = eta(2);
tauL = eta(3);
azR = eta(4);
elR = eta(5);
tauR = eta(6);
vartheta2 = eta(7);
vartheta3 = eta(8);

% setup
if ~exist('oNum', 'var') || ~exist('dNum', 'var')
    oNum = 100;
    dNum = 100;
end

% determine a distance range of AOD line
AOD_L = [azL, elL].';
tAOD_L = get_dir_from_angle(AOD_L);
tAOD_L_global = Rb*tAOD_L;
AOD_R = [azR, elR].';
tAOD_R = get_dir_from_angle(AOD_R);
tAOD_R_global = Rb*tAOD_R;
o_range = linspace(olim(1), olim(2), oNum);
o_range = o_range + randn(1)*(o_range(2)-o_range(1));   % add a small random offset
d_range = linspace(dlim(1), dlim(2), dNum );
d_range = d_range + randn(1)*(d_range(2)-d_range(1));   % add a small random offset
cost = zeros(oNum, dNum);

for i = 1:oNum
    o3 = o_range(i);
    for j = 1:dNum
        clkb = d_range(j);
        d0 = tauL - clkb;
        pu = pb + tAOD_L_global*d0;
        dr = tauR - clkb;
        pr = get_intersection_ellipsoid_line(pb,pu,dr,tAOD_R_global);
        
        [varth2, varth3] = get_vartheta(o3,pu,pr,pb,sp);
        cost(i,j) = sqrt( (vartheta2-varth2)^2 + (vartheta3-varth3)^2 );
    end
end
[~, ind_o] = min(cost);
[~, ind_d] = min(min(cost));

clkb_est = d_range(ind_d);
o3_est = o_range(ind_o(ind_d));


% generate localization parameters
d0 = tauL - clkb_est;
pu = pb + tAOD_L_global*d0;
dr = tauR - clkb_est;
pr = get_intersection_ellipsoid_line(pb,pu,dr,tAOD_R_global);
r_est = [pu; pr; o3_est; clkb_est];



%% ================= functions
    function pr = get_intersection_ellipsoid_line(pb,pu,dr,tg)
        x = ( dr^2 - norm(pb-pu,2)^2 ) / ( 2*(dr + tg.'*(pb-pu) ) );
        pr = pb + x*tg; 
    end

    function [vartheta2, vartheta3] = get_vartheta(o3,Pu,Pr,Pb,sp)
        dbr = norm(Pb-Pr,2);
        dur = norm(Pu-Pr,2);
        % intermediate unknowns
        Or_Euler = [o3;sp.Or_Euler(2:3)];
        Rr = eul2rotm(reshape(deg2rad(Or_Euler), [1,3]), 'ZYX');
        tRU = Rr.'*( (Pu - Pr)/dur);
        tRB = Rr.'*( (Pb - Pr)/dbr);
        vartheta2 = tRU(2) + tRB(2);
        vartheta3 = tRU(3) + tRB(3);
    end
end


