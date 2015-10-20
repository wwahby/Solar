%% Calculate excess hole concentration in silicon under solar illumination

lambda_min_nm = 1; % (nm) Minimum wavelength to consider
lambda_max_nm = 1500; % (nm) Maximum wavelength to consider
Npts_lambda = 16; % Number of points for wavelength
Npts_x = 1e2; % Number of points for x dimension
Npts_z = 1e2; % Number of points for z dimension
m_max = 10; % number of terms to use in fourier expansions
W_m = 20e-6; % (m) Horizontal separation between n-doped regions
d_m = 1000e-6; % (m) Depth of cell
T = 300; % (K)
tau_n = 1e-6; % (s)
Ln_m = 33e-6; % (m)
np = 1e4; % (cm^-3);


% Convert all units to cm
Ln_cm = Ln_m * 1e2;
W_cm = W_m  *1e2;
d_cm = d_m * 1e2;
Dn_cm2 = Ln_cm^2/tau_n;
q = 1.602e-19; % (C) Electron charge
Jo = q*Dn_cm2/Ln_cm*np;

zvec_m = logspace(-9,-3,Npts_z);
zvec_cm = zvec_m * 1e2;

flux_mat = zeros(Npts_z, Npts_lambda);

for zind = 1:length(zvec_cm)
    z_cm = zvec_cm(zind);
    lambda_vec_nm = linspace(lambda_min_nm, lambda_max_nm, Npts_lambda);
    [eta_col_vec, eta_abs_vec, flux_vec_m2, alpha_vec_cm, flux_wehrli_vec_m2] = vsc.calc_collection_efficiency(lambda_vec_nm, m_max, W_cm, d_cm, Dn_cm2, Ln_cm, Sn_cm);
    flux_vec_cm2 = flux_vec_m2/1e4;
    flux_wehrli_vec_cm2 = flux_wehrli_vec_m2/1e4;
    
    flux_mat(zind, :) = flux_vec_cm2 .* exp(-alpha_vec_cm .*z_cm);
end


%% Plots
colors = {'k', 'b', 'r', 'g'};
linestyles = {'-', '--', '-.', ':'};

f1 = figure(1);
clf
hold on

f2 = figure(2);
clf
hold on

z_diff_vec = diff(zvec_m);
for lind = 1:length(lambda_vec_nm)
    col_ind = mod(lind-1, length(colors))+1;
    linestyle_ind = floor( (lind-1)/(length(colors)) ) + 1;
    flux_z_vec = zeros(1, Npts_z);
    flux_z_vec(:) = flux_mat(:, lind)';
    flux_z_diff_vec = diff(flux_z_vec * 1e4);
    flux_z_deriv_vec = flux_z_diff_vec./z_diff_vec;
    figure(1)
    plot( zvec_m*1e6, flux_z_vec, 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
    figure(2)
    plot( z_diff_vec*1e6, -flux_z_deriv_vec, 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
end

figure(1)
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e0 1e20])
xlabel('Depth (um)')
ylabel('Photon Flux (cm^{-2})')
fixfigs(1,3,14,12)

figure(2)
set(gca,'yscale','log')
set(gca,'xscale','log')
%ylim([1e0 1e20])
xlabel('Depth (um)')
ylabel('Derivative of Photon Flux (m^{-1})')
fixfigs(2,3,14,12)
