lambda_nm_vec = linspace(1,1100, 1e2);


r_dep_m = 2e-6;
L_col_m = 10e-6;
r_col_m = r_dep_m + L_col_m;
zo_m = 12e-6;

G_dep_per_nm_vec = zeros(1,length(lambda_nm_vec));
G_col_per_nm_vec = zeros(1,length(lambda_nm_vec));
flux_per_m2nm_vec = zeros(1, length(lambda_nm_vec));
irradiance_per_m2nm_vec = zeros(1, length(lambda_nm_vec));
for lind = 1:length(lambda_nm_vec)
    lambda_nm = lambda_nm_vec(lind);
    [G_dep_per_nm, flux_per_m2nm, irradiance_per_m2nm, alpha_cm] = cpm.calc_generation_rate_in_sphere(lambda_nm, r_dep_m, zo_m);
    [G_col_per_nm, flux_per_m2nm, irradiance_per_m2nm, alpha_cm] = cpm.calc_generation_rate_in_sphere(lambda_nm, r_col_m, zo_m);
    
    G_dep_per_nm_vec(lind) = G_dep_per_nm;
    G_col_per_nm_vec(lind) = G_col_per_nm;
    
    flux_per_m2nm_vec(lind) = flux_per_m2nm;
    irradiance_per_m2nm_vec(lind) = irradiance_per_m2nm;
end


q = 1.602e-19; % (C) Electron charge
E_ev = 1240./lambda_nm_vec;
E_J = E_ev*q;

G_dep_tot = trapz(lambda_nm_vec, G_dep_per_nm_vec);
G_col_tot = trapz(lambda_nm_vec, G_col_per_nm_vec);
flux_tot_per_m2 = trapz(lambda_nm_vec, flux_per_m2nm_vec);

G_act = G_col_tot - G_dep_tot;


L_side = 2*r_col_m;
A = L_side^2;

i_act = G_act * q;
J_act = i_act/A/1e4; % current per surface area (cm2)

flux_tot = flux_tot_per_m2 * A;
eff_col = G_act/flux_tot

%%
figure(3)
clf
hold on
semilogy(lambda_nm_vec, G_dep_per_nm_vec, 'b')
semilogy(lambda_nm_vec, G_col_per_nm_vec, 'r')
ylim( [1e-10 10*max( [max(G_dep_per_nm_vec) max(G_col_per_nm_vec)]) ])
set(gca,'yscale', 'log')
fixfigs(3,3,14,12)

figure(4)
clf
hold on
plot(lambda_nm_vec, G_dep_per_nm_vec, 'b')
plot(lambda_nm_vec, G_col_per_nm_vec, 'r')
fixfigs(4,3,14,12)

figure(5)
clf
hold on
plot(lambda_nm_vec, flux_per_m2nm_vec,'b')
set(gca,'yscale','log')
xlabel('Wavelength (nm)')
ylabel('Spectral Flux Density (m^{-2}nm^{-1})')
fixfigs(5,3,14,12)

figure(6)
clf
hold on
plot(lambda_nm_vec, irradiance_per_m2nm_vec, 'b')
set(gca,'yscale','log')
xlabel('Wavelength (nm)')
ylabel('Spectral Irradiance (Wm^{-2}nm^{-1})')
fixfigs(6,3,14,12)