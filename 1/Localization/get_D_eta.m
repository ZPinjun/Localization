function D_eta = get_D_eta(q,x,hmn_a,aBm,aSn,fk,fc,c,eta,pAE_B,pAE_S)

theta_az = eta(1);
theta_el = eta(2);
phi_az = eta(3);
phi_el = eta(4);
tau = eta(5)/c;
hmn_p = 2*pi*fc*tau;
hmn = hmn_a * exp(-1j*hmn_p);
df = fk-fc;
zeta_n = hmn*exp(-1j*2*pi*df*tau)*q*aSn*aBm.'*x;


t_AOD_tilde = [ -cosd(theta_el)*sind(theta_az);
                cosd(theta_el)*cosd(theta_az);
                0 ];
t_AOD_check = [ -sind(theta_el)*cosd(theta_az);
                -sind(theta_el)*sind(theta_az);
                cosd(theta_el) ];
t_AOA_tilde = [ -cosd(phi_el)*sind(phi_az);
                cosd(phi_el)*cosd(phi_az);
                0 ];
t_AOA_check = [ -sind(phi_el)*cosd(phi_az);
                -sind(phi_el)*sind(phi_az);
                cosd(phi_el) ];


aBm_tilde = 1j*2*pi*fk/c * aBm .* (pAE_B.'*t_AOD_tilde);
aBm_check = 1j*2*pi*fk/c * aBm .* (pAE_B.'*t_AOD_check);
aSn_tilde = 1j*2*pi*fk/c * aSn .* (pAE_S.'*t_AOA_tilde);
aSn_check = 1j*2*pi*fk/c * aSn .* (pAE_S.'*t_AOA_check);


D_theta_az = hmn * exp(-1j*2*pi*df*tau) * q*aSn*aBm_tilde.'*x;
D_theta_el = hmn * exp(-1j*2*pi*df*tau) * q*aSn*aBm_check.'*x;
D_phi_az = hmn * exp(-1j*2*pi*df*tau) * q*aSn_tilde*aBm.'*x;
D_phi_el = hmn * exp(-1j*2*pi*df*tau) * q*aSn_check*aBm.'*x;
D_tau = -1j*2*pi*df/c * zeta_n;   % transfer tau in [m]
D_hmn_a = zeta_n / hmn_a;
D_hmn_p = -1j * zeta_n;

D_eta = [D_theta_az, D_theta_el, D_phi_az, D_phi_el, D_tau, D_hmn_a, D_hmn_p];

end

