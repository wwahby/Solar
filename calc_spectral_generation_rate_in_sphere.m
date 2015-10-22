function [G_per_nm, flux_per_m2nm, irradiance_per_m2nm, alpha_cm] = calc_spectral_generation_rate_in_sphere(lambda_nm, r_m, zo_m)
%% Calculates carrier generation rate in a sphere of radius r, centered on depth zo, for illumination of wavelength lambda


[flux_per_m2, flux_wehrli_per_m2, flux_per_m2nm, flux_wehrli_per_m2nm, irradiance_per_m2nm, irradiance_wehrli_per_m2nm] = vsc.get_photon_flux_per_m2( lambda_nm );
alpha_cm = vsc.get_absorption_coefficient(lambda_nm);

%% Conversions
r = r_m*1e2;
zo = zo_m*1e2;
alpha = alpha_cm;

za = zo-r;
zb = zo+r;

phi_o_cm2nm = flux_per_m2nm/1e4;
phi_o = phi_o_cm2nm;

%% Calculations
% Simplified expression for full sphere -- not necessary since partial
% sphere expression is also very compact
%G_per_nm  = 2*pi/alpha^2*phi_o * ( exp(-alpha*zb)*(alpha*r + 1) + exp(-alpha*za)*(alpha*r-1) ); % total generation rate per wavelength in spherical vol

% Clip collection sphere if it extends past surface
if (zo < r)
    za = 0;
end

G_per_nm_A = exp(-alpha*za)*(alpha^2*(r^2-(zo-za)^2) + 2*alpha*(zo-za)-2);
G_per_nm_B = exp(-alpha*zb)*(alpha^2*(r^2-(zo-zb)^2) + 2*alpha*(zb-zo)+2);
G_per_nm = pi/alpha^2*phi_o * (G_per_nm_A  + G_per_nm_B);