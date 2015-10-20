%function nvol = calc_nvol(phi_o_cm2, r_m, alpha_cm, zo_m)

lambda_nm = 500;
r_m = 2e-6;
zo_m = 20e-6;

[flux_per_m2, flux_wehrli_per_m2] = vsc.get_photon_flux_per_m2( lambda_nm );
alpha_cm = vsc.get_absorption_coefficient(lambda_nm);

phi_o_cm2 = flux_per_m2 /1e4;

%% Conversions
r = r_m*1e2;
zo = zo_m*1e2;
alpha = alpha_cm;
phi_o = phi_o_cm2;

zt = zo-r;
zb = zo+r;

%% Calculations
nvol_a = pi*r^2*phi_o*(exp(-alpha*zt) - exp(-alpha*zb));
nvol_b = pi/alpha^2*phi_o*exp(-alpha*zo);
nvol_c = exp(-alpha*r)*(alpha^2*r^2 + 2*alpha*r + 2);
nvol_d = exp(alpha*r) *(alpha^2*r^2 - 2*alpha*r + 2);

%% Combine
nvol = nvol_a + nvol_b*(nvol_c - nvol_d)