lambda_min = 400; % (nm) Minimum wavelength to consider
lambda_max = 1100; % (nm) Maximum wavelength to consider
Npts_lambda = 1e1; % Number of points for wavelength
Npts_x = 1e2; % Number of points for x dimension
Npts_z = 1e2; % Number of points for z dimension
m_max = 50; % number of terms to use in fourier expansions
W_m = 20e-6; % (m) Horizontal separation between n-doped regions
d_m = 1000e-6; % (m) Depth of cell
T = 300; % (K)
tau_n = 1e6; % (s)
Ln_m = 33e-6; % (m)
Jo = 1e-6; % ??

Sn = 3000;
Sn_top = Sn;
Sn_bot = Sn;

Ln = Ln_m * 1e2;
W = W_m * 1e2;
d = d_m * 1e2;
Dn = Ln^2/tau_n;

%     Isc = calc_Isc(m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, Ln, Sn);
%     Voc_z = calc_Voc_z(z, m_max, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot, T, Jo);

lambda_vec = linspace(lambda_min, lambda_max, Npts_lambda);
[eta_col_vec, eta_abs_vec, flux_vec, alpha_vec] = vsc.calc_collection_efficiency(lambda_vec, m_max, W, d, Dn, Ln, Sn);



d_m_vec = linspace(1e-6, 1000e-6, 1e3);
d_vec = d_m_vec * 1e2; % (cm)

Voc_vec = zeros(1, length(d_vec));
Isc_vec = zeros(1, length(d_vec));
for dind = 1:length(d_vec)
    d = d_vec(dind);
    Voc = vsc.calc_Voc_z(d, m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot, T, Jo);
    Isc = vsc.calc_Isc(m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, Ln, Sn);
    Voc_vec(dind) = Voc;
    Isc_vec(dind) = Isc;
end

figure(1)
clf
hold on
plot(lambda_vec, eta_col_vec,'b')
plot(lambda_vec, eta_abs_vec,'r')
set(gca, 'yscale','log')
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
