function [G_shell, G_inner, G_outer, flux_tot, eff_col_shell, flux_tot_per_m2, i_act, J_act] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_m)

G_inner_per_nm_vec = zeros(1,length(lambda_nm_vec));
G_outer_per_nm_vec = zeros(1,length(lambda_nm_vec));
flux_per_m2nm_vec = zeros(1, length(lambda_nm_vec));
irradiance_per_m2nm_vec = zeros(1, length(lambda_nm_vec));
for lind = 1:length(lambda_nm_vec)
    lambda_nm = lambda_nm_vec(lind);
    [G_inner_per_nm, flux_per_m2nm, irradiance_per_m2nm, alpha_cm] = calc_generation_rate_in_sphere(lambda_nm, r_inner_m, zo_m);
    [G_outer_per_nm, flux_per_m2nm, irradiance_per_m2nm, alpha_cm] = calc_generation_rate_in_sphere(lambda_nm, r_outer_m, zo_m);
    
    G_inner_per_nm_vec(lind) = G_inner_per_nm;
    G_outer_per_nm_vec(lind) = G_outer_per_nm;
    
    flux_per_m2nm_vec(lind) = flux_per_m2nm;
    irradiance_per_m2nm_vec(lind) = irradiance_per_m2nm;
end


q = 1.602e-19; % (C) Electron charge
E_ev = 1240./lambda_nm_vec;
E_J = E_ev*q;

G_inner = trapz(lambda_nm_vec, G_inner_per_nm_vec);
G_outer = trapz(lambda_nm_vec, G_outer_per_nm_vec);
flux_tot_per_m2 = trapz(lambda_nm_vec, flux_per_m2nm_vec);

G_shell = G_outer - G_inner;


L_side = 2*r_outer_m;
A = L_side^2;

i_act = G_shell * q;
J_act = i_act/A/1e4; % current per surface area (cm2)

flux_tot = flux_tot_per_m2 * A;
eff_col_shell = G_shell/flux_tot;