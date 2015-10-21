lambda_nm_vec = linspace(1,1100, 1e2);


r_dep_m = 1e-6;
L_col_m = 2e-6;
r_col_m = r_dep_m + L_col_m;
zo_m = 10e-6;

nvol_dep_vec = zeros(1,length(lambda_nm_vec));
nvol_dep_per_nm_vec = zeros(1,length(lambda_nm_vec));
G_dep_per_nm_vec = zeros(1,length(lambda_nm_vec));
nvol_col_vec = zeros(1,length(lambda_nm_vec));
nvol_col_per_nm_vec = zeros(1,length(lambda_nm_vec));
G_col_per_nm_vec = zeros(1,length(lambda_nm_vec));
flux_per_m2nm_vec = zeros(1, length(lambda_nm_vec));
irradiance_per_m2nm_vec = zeros(1, length(lambda_nm_vec));
for lind = 1:length(lambda_nm_vec)
    lambda_nm = lambda_nm_vec(lind);
    [nvol_dep, nvol_dep_per_nm, G_dep_per_nm, flux_per_m2, flux_per_m2nm, irradiance_per_m2nm] = calc_nvol(lambda_nm, r_dep_m, zo_m);
    [nvol_col, nvol_col_per_nm, G_col_per_nm, flux_per_m2, flux_per_m2nm, irradiance_per_m2nm] = calc_nvol(lambda_nm, r_col_m, zo_m);
    nvol_dep_vec(lind) = nvol_dep;
    nvol_dep_per_nm_vec(lind) = nvol_dep_per_nm;
    G_dep_per_nm_vec(lind) = G_dep_per_nm;
    nvol_col_vec(lind) = nvol_col;
    nvol_col_per_nm_vec(lind) = nvol_col_per_nm;
    G_col_per_nm_vec(lind) = G_col_per_nm;
    
    flux_per_m2nm_vec(lind) = flux_per_m2nm;
    irradiance_per_m2nm_vec(lind) = irradiance_per_m2nm;
end


q = 1.602e-19; % (C) Electron charge
E_ev = 1240./lambda_nm_vec;
E_J = E_ev*q;
    
    
nvol_dep_tot = trapz(lambda_nm_vec, nvol_dep_vec);
nvol_dep_tot_2 = trapz(lambda_nm_vec, nvol_dep_per_nm_vec);
nvol_dep_tot_3 = trapz(lambda_nm_vec, nvol_dep_vec./lambda_nm_vec);
G_dep_tot = trapz(lambda_nm_vec, G_dep_per_nm_vec)

nvol_col_tot = trapz(lambda_nm_vec, nvol_col_vec);
nvol_col_tot_2 = trapz(lambda_nm_vec, nvol_col_per_nm_vec);
nvol_col_tot_3 = trapz(lambda_nm_vec, nvol_col_vec./lambda_nm_vec);
G_col_tot = trapz(lambda_nm_vec, G_col_per_nm_vec)

flux_tot_per_m2 = trapz(lambda_nm_vec, flux_per_m2nm_vec)

nvol_act = nvol_col_tot_2 - nvol_dep_tot_2
G_act = G_col_tot - G_dep_tot

i_act = nvol_act * q;
L_side = 2*r_col_m;
A = L_side^2;
flux_tot = flux_tot_per_m2 * A;
eff_col = nvol_act/flux_tot
eff_col_2 = G_act/flux_tot
J_act = i_act/A/1e4;

%%
figure(1)
clf
hold on
semilogy(lambda_nm_vec, nvol_dep_vec, 'b')
semilogy(lambda_nm_vec, nvol_col_vec, 'r')
ylim( [1e-10 10*max( [max(nvol_dep_vec) max(nvol_col_vec)]) ])
fixfigs(1,3,14,12)

figure(2)
clf
hold
plot(lambda_nm_vec, nvol_dep_vec, 'b')
plot(lambda_nm_vec, nvol_col_vec, 'r')
fixfigs(2,3,14,12)

figure(3)
clf
hold on
semilogy(lambda_nm_vec, nvol_dep_per_nm_vec, 'b')
semilogy(lambda_nm_vec, nvol_col_per_nm_vec, 'r')
ylim( [1e-10 10*max( [max(nvol_dep_per_nm_vec) max(nvol_col_per_nm_vec)]) ])
fixfigs(3,3,14,12)

figure(4)
clf
hold on
plot(lambda_nm_vec, nvol_dep_per_nm_vec, 'b')
plot(lambda_nm_vec, nvol_col_per_nm_vec, 'r')
fixfigs(4,3,14,12)

figure(5)
clf
hold on
plot(lambda_nm_vec, flux_per_m2nm_vec,'b')
set(gca,'yscale','log')
xlabel('Wavelength (nm)')
ylabel('Spectral Flux (m^{-2}nm^{-1})')
fixfigs(5,3,14,12)

figure(6)
clf
hold on
plot(lambda_nm_vec, irradiance_per_m2nm_vec, 'b')
set(gca,'yscale','log')
xlabel('Wavelength (nm)')
ylabel('Spectral Irradiance (Wm^{-2}nm^{-1})')
fixfigs(6,3,14,12)