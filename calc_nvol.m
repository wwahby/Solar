function [nvol, nvol_per_nm] = calc_nvol(lambda_nm, r_m, zo_m)


[flux_per_m2, flux_wehrli_per_m2, irradiance_per_m2nm, irradiance_wehrli_per_m2nm] = vsc.get_photon_flux_per_m2( lambda_nm );
alpha_cm = vsc.get_absorption_coefficient(lambda_nm);

phi_o_cm2 = flux_per_m2/1e4;

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

%% Calculations
phi_o_cm2nm = irradiance_per_m2nm/1e4;
phi_o = phi_o_cm2nm;

nvol_per_nm_a = pi*r^2*phi_o*(exp(-alpha*zt) - exp(-alpha*zb));
nvol_per_nm_b = pi/alpha^2*phi_o*exp(-alpha*zo);
nvol_per_nm_c = exp(-alpha*r)*(alpha^2*r^2 + 2*alpha*r + 2);
nvol_per_nm_d = exp(alpha*r) *(alpha^2*r^2 - 2*alpha*r + 2);

%% Combine
nvol = nvol_a + nvol_b*(nvol_c - nvol_d);
nvol_per_nm_ = nvol_per_nm_a + nvol_per_nm_b*(nvol_per_nm_c - nvol_per_nm_d);