lambda_min_nm = 400; % (nm) Minimum wavelength to consider
lambda_max_nm = 1100; % (nm) Maximum wavelength to consider
Npts_lambda = 1e2; % Number of points for wavelength
Npts_x = 1e2; % Number of points for x dimension
Npts_z = 1e2; % Number of points for z dimension
m_max = 10; % number of terms to use in fourier expansions
W_m = 20e-6; % (m) Horizontal separation between n-doped regions
d_m = 1000e-6; % (m) Depth of cell
T = 300; % (K)
tau_n = 1e-6; % (s)
Ln_m = 33e-6; % (m)
np = 1e4; % (cm^-3);

sn_vec_cm = [0, 1000, 2000, 3000];

% Convert all units to cm
Ln_cm = Ln_m * 1e2;
W_cm = W_m  *1e2;
d_cm = d_m * 1e2;
Dn_cm2 = Ln_cm^2/tau_n;
q = 1.602e-19; % (C) Electron charge
Jo = q*Dn_cm2/Ln_cm*np;


%% D75 Fig2
figure(1)
clf
hold on
colors = {'k', 'b', 'r', 'g', 'm', 'c'};
linestyles = {'-', '--', '-.', ':'};
for snind = 1:length(sn_vec_cm)
    Sn_cm = sn_vec_cm(snind);
    Sn_top_cm = Sn_cm;
    Sn_bot_cm = Sn_cm;
    
    lambda_vec_nm = linspace(lambda_min_nm, lambda_max_nm, Npts_lambda);
    [eta_col_vec, eta_abs_vec, flux_vec_m2, alpha_vec_cm] = vsc.calc_collection_efficiency(lambda_vec_nm, m_max, W_cm, d_cm, Dn_cm2, Ln_cm, Sn_cm);
    
    plot(lambda_vec_nm, eta_col_vec, 'color', 'b', 'linestyle', linestyles{snind} )
    plot(lambda_vec_nm, eta_abs_vec, 'color', 'r', 'linestyle', linestyles{snind} )
end



%set(gca, 'yscale','log')
ylim([0.80 1.0])
xlabel('Wavelength (nm)')
ylabel('Efficiency')
title('Dhariwal 1975 - Figure 2')
fixfigs(1,3,14,12)

figure(2)
clf
plot(lambda_vec_nm, flux_vec_m2)
xlabel('Wavelength (nm)')
ylabel('Photon flux')
fixfigs(2,3,14,12)

figure(3)
clf
plot(lambda_vec_nm, alpha_vec_cm)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (1/cm)')
set(gca, 'yscale','log')
fixfigs(3,3,14,12)

%% D75 Figure 3
d_m = 1000e-6;
d_cm = d_m*100;
Ln_m = 33e-6;
Ln_cm = Ln_m*100;
Sn_cm = 0;
tau_n = 1e-6; % (s)
Dn_cm2 = Ln_cm^2/tau_n;
W_cm = 1;

zvec_m = logspace(-9,-3,Npts_z);
%zvec_m = linspace(0, 1000e-6, Npts_z);
zvec_cm = zvec_m * 1e2;
gen_vec_z = zeros(1, length(zvec_cm));
gen_wehrli_vec_z = zeros(1, length(zvec_cm));
for zind = 1:length(zvec_cm)
    z_cm = zvec_cm(zind);
    lambda_vec_nm = linspace(0, lambda_max_nm, Npts_lambda);
    [eta_col_vec, eta_abs_vec, flux_vec_m2, alpha_vec_cm, flux_wehrli_vec_m2] = vsc.calc_collection_efficiency(lambda_vec_nm, m_max, W_cm, d_cm, Dn_cm2, Ln_cm, Sn_cm);
    flux_vec_cm2 = flux_vec_m2/1e4;
    flux_wehrli_vec_cm2 = flux_wehrli_vec_m2/1e4;
    
    gen_vec_integrand = flux_vec_cm2 .* alpha_vec_cm .* exp(-alpha_vec_cm * z_cm);
    gen_wehrli_vec_integrand = flux_wehrli_vec_cm2 .* alpha_vec_cm .* exp(-alpha_vec_cm * z_cm);
    lambda_vec_cm = lambda_vec_nm /1e7;
    gen_vec_z(zind) = trapz(lambda_vec_cm, gen_vec_integrand);
    gen_wehrli_vec_z(zind) = trapz(lambda_vec_cm, gen_wehrli_vec_integrand);
end

figure(4)
clf
hold on
plot(zvec_m*1e6, gen_vec_z)
plot(zvec_m*1e6, gen_wehrli_vec_z,'r--')
set(gca,'yscale','log')
grid on
%ylim([1e16 1e24])
xlabel('Vertical Distance (um)')
ylabel('Generation Constant (cm^{-3}s^{-1})')
fixfigs(4,3,14,12)

    


%%
d_max_m = 1000e-6;
d_m_vec = linspace(1e-6, d_max_m, 1e2);
d_vec = d_m_vec * 1e2; % (cm)
d_m_fixed = d_max_m;
d_fixed = d_m_fixed * 1e2;

Voc_vec = zeros(1, length(d_vec));
Jsc_vec = zeros(1, length(d_vec));
Isc_vec = zeros(1, length(d_vec));
Sn_top_cm = 0;
Sn_bot_cm = 0;
for dind = 1:length(d_vec)
    d_cm = d_vec(dind);
    Jsc_z = vsc.calc_Jsc_z(d_cm, m_max, lambda_min_nm, lambda_max_nm, Npts_lambda, W_cm, d_cm, Dn_cm2, tau_n, Sn_top_cm, Sn_bot_cm);
    Voc = vsc.calc_Voc_z(d_cm, m_max, lambda_min_nm, lambda_max_nm, Npts_lambda, W_cm, d_cm, Dn_cm2, tau_n, Sn_top_cm, Sn_bot_cm, T, Jo);
    Isc = vsc.calc_Isc(m_max, lambda_min_nm, lambda_max_nm, Npts_lambda, W_cm, d_cm, Dn_cm2, Ln_cm, Sn_cm);
    Voc_vec(dind) = Voc;
    Isc_vec(dind) = Isc;
    Jsc_vec(dind) = Jsc_z;
end

%% Power out vs power in
photon_flux_cm2_vec = zeros(1,length(lambda_vec_nm));
photon_power_cm2_vec = zeros(1, length(lambda_vec_nm));
photon_irradiance_cm2nm_vec = zeros(1, length(lambda_vec_nm));
for lind = 1:length(lambda_vec_nm)
    lambda = lambda_vec_nm(lind);
    flux_per_m2 = vsc.get_photon_flux_per_m2( lambda );
    photon_flux_cm2 = flux_per_m2/100^2; % N is per cm^2
    photon_flux_cm2_vec(lind) = photon_flux_cm2;
    [photon_power_m2, photon_irradiance_m2nm] = vsc.get_photon_power_per_m2( lambda );
    photon_power_cm2_vec(lind) = photon_power_m2 / 100^2;
    photon_irradiance_cm2nm_vec(lind) = photon_irradiance_m2nm / 100^2;
end

E_photon_ev = 1240./lambda_vec_nm;
E_photon = E_photon_ev*q;

Pout = Voc_vec .* Isc_vec; % Power output per cm of solar cell length (function of wafer thickness
Pin = sum(photon_power_cm2_vec) .* d_vec; % Input power per cm of solar cell length
%Pin = trapz(lambda_vec, photon_irradiance_cm2nm_vec) * d_vec;

eta_eff = Pout./Pin;

%%
figure(5)
clf
plot(d_m_vec*1e6, Isc_vec)
xlabel('Thickness (microns)')
ylabel('Isc (A)')
%set(gca, 'yscale','log')
fixfigs(5,3,14,12)

figure(6)
clf
plot(d_m_vec*1e6, Voc_vec)
xlabel('Thickness (microns)')
ylabel('Voc (V)')
%set(gca, 'yscale','log')
fixfigs(6,3,14,12)

figure(7)
clf
hold on
plot(d_m_vec*1e6, eta_eff)
xlabel('Thickness (microns)')
ylabel('Efficiency')
%ylim([0 1])
fixfigs(7,3,14,12)
