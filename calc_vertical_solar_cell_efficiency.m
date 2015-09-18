lambda_min = 400; % (nm) Minimum wavelength to consider
lambda_max = 1100; % (nm) Maximum wavelength to consider
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

sn_vec = [0, 1000, 2000, 3000];

% Convert all units to cm
Ln = Ln_m * 1e2;
W = W_m  *1e2;
d = d_m * 1e2;
Dn = Ln^2/tau_n;
q = 1.602e-19; % (C) Electron charge
Jo = q*Dn/Ln*np;



figure(1)
clf
hold on
colors = {'k', 'b', 'r', 'g', 'm', 'c'};
linestyles = {'-', '--', '-.', ':'};
for snind = 1:length(sn_vec)
    Sn = sn_vec(snind);
    Sn_top = Sn;
    Sn_bot = Sn;
    
    lambda_vec = linspace(lambda_min, lambda_max, Npts_lambda);
    [eta_col_vec, eta_abs_vec, flux_vec, alpha_vec] = vsc.calc_collection_efficiency(lambda_vec, m_max, W, d, Dn, Ln, Sn);
    plot(lambda_vec, eta_col_vec, 'color', 'b', 'linestyle', linestyles{snind} )
    plot(lambda_vec, eta_abs_vec, 'color', 'r', 'linestyle', linestyles{snind} )
end



%set(gca, 'yscale','log')
ylim([0.80 1.0])
xlabel('Wavelength (nm)')
ylabel('Efficiency')
fixfigs(1,3,14,12)

figure(2)
clf
plot(lambda_vec, flux_vec)
xlabel('Wavelength (nm)')
ylabel('Photon flux')
fixfigs(2,3,14,12)

figure(3)
clf
plot(lambda_vec, alpha_vec)
xlabel('Wavelength (nm)')
ylabel('Absorption Coefficient (1/cm)')
set(gca, 'yscale','log')
fixfigs(3,3,14,12)


%%
d_max_m = 1000e-6;
d_m_vec = linspace(1e-6, d_max_m, 1e2);
d_vec = d_m_vec * 1e2; % (cm)
d_m_fixed = d_max_m;
d_fixed = d_m_fixed * 1e2;

Voc_vec = zeros(1, length(d_vec));
Jsc_vec = zeros(1, length(d_vec));
Isc_vec = zeros(1, length(d_vec));
Sn_top = 0;
Sn_bot = 0;
for dind = 1:length(d_vec)
    d = d_vec(dind);
    Jsc_z = vsc.calc_Jsc_z(d, m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot);
    Voc = vsc.calc_Voc_z(d, m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot, T, Jo);
    Isc = vsc.calc_Isc(m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, Ln, Sn);
    Voc_vec(dind) = Voc;
    Isc_vec(dind) = Isc;
    Jsc_vec(dind) = Jsc_z;
end

%% Power out vs power in
photon_flux_cm2_vec = zeros(1,length(lambda_vec));
photon_power_cm2_vec = zeros(1, length(lambda_vec));
photon_irradiance_cm2nm_vec = zeros(1, length(lambda_vec));
for lind = 1:length(lambda_vec)
    lambda = lambda_vec(lind);
    flux_per_m2 = vsc.get_photon_flux_per_m2( lambda );
    photon_flux_cm2 = flux_per_m2/100^2; % N is per cm^2
    photon_flux_cm2_vec(lind) = photon_flux_cm2;
    [photon_power_m2, photon_irradiance_m2nm] = vsc.get_photon_power_per_m2( lambda );
    photon_power_cm2_vec(lind) = photon_power_m2 / 100^2;
    photon_irradiance_cm2nm_vec(lind) = photon_irradiance_m2nm / 100^2;
end

E_photon_ev = 1240./lambda_vec;
E_photon = E_photon_ev*q;

Pout = Voc_vec .* Isc_vec; % Power output per cm of solar cell length (function of wafer thickness
Pin = sum(photon_power_cm2_vec) .* d_vec; % Input power per cm of solar cell length
%Pin = trapz(lambda_vec, photon_irradiance_cm2nm_vec) * d_vec;

eta_eff = Pout./Pin;

%%
figure(4)
clf
plot(d_m_vec*1e6, Isc_vec)
xlabel('Thickness (microns)')
ylabel('Isc (A)')
%set(gca, 'yscale','log')
fixfigs(4,3,14,12)

figure(5)
clf
plot(d_m_vec*1e6, Voc_vec)
xlabel('Thickness (microns)')
ylabel('Voc (V)')
%set(gca, 'yscale','log')
fixfigs(5,3,14,12)

figure(6)
clf
hold on
plot(d_m_vec*1e6, eta_eff)
xlabel('Thickness (microns)')
ylabel('Efficiency')
%ylim([0 1])
fixfigs(6,3,14,12)
