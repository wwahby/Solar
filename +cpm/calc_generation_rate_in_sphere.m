function [G_per_nm, flux_per_m2nm, irradiance_per_m2nm, alpha_cm] = calc_generation_rate_in_sphere(lambda_nm, r_m, zo_m)
%% Calculates carrier generation rate in a sphere of radius r, centered on depth zo, for illumination of wavelength lambda
% Assumes zo-r > 0

[flux_per_m2, flux_wehrli_per_m2, flux_per_m2nm, flux_wehrli_per_m2nm, irradiance_per_m2nm, irradiance_wehrli_per_m2nm] = vsc.get_photon_flux_per_m2( lambda_nm );
alpha_cm = vsc.get_absorption_coefficient(lambda_nm);

%% Conversions
r = r_m*1e2;
zo = zo_m*1e2;
alpha = alpha_cm;

za = zo-r;
zb = zo+r;

%% Calculations
phi_o_cm2nm = flux_per_m2nm/1e4;
phi_o = phi_o_cm2nm;
G_per_nm  = 2*pi/alpha^2*phi_o * ( exp(-alpha*zb)*(alpha*r + 1) + exp(-alpha*za)*(alpha*r-1) ); % total generation rate per wavelength in spherical vol
