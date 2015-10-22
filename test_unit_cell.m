lambda_nm_vec = linspace(1,1100, 1e2);

r_dep_m = 1e-6;
L_col_m = 2e-6;
r_col_m = r_dep_m + L_col_m;
zo_m = 4e-6;

r_inner_m = r_dep_m;
r_outer_m = r_col_m;
zo_top_m = zo_m;
zo_bot_m = zo_top_m + sqrt(3)*r_outer_m;
[G_shell_top, G_inner_top, G_outer_top, flux_tot, eff_col_shell_top, flux_tot_per_m2, i_act_top, J_act_top] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_top_m);
[G_shell_bot, G_inner_bot, G_outer_bot, flux_tot, eff_col_shell_bot, flux_tot_per_m2, i_act_bot, J_act_bot] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_bot_m);

G_tot = G_shell_top + G_shell_bot;
i_act_tot = i_act_top + i_act_bot;
J_act_tot = J_act_top + J_act_bot;



num_collectors = 10;


zo_m_vec = (0:4) * 1e-6;
figure(1)
clf
hold on
colors = {'k', 'b', 'g', 'y', 'r' };
for zind = 1:length(zo_m_vec)
    Gtot_vec = zeros(1,num_collectors);
    zo_m = zo_m_vec(zind);
    for nind = 1:num_collectors
        depth_m = zo_m  + (nind)*r_outer_m*sqrt(2);
        [G_shell, G_inner, G_outer, flux_tot, eff_col_shell, flux_tot_per_m2, i_act, J_act] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, depth_m);
        
        Gtot_vec(nind) = 2*G_shell; % fcc lattice -- two shells per layer

    end
    Gtot = sum(Gtot_vec);
    
    L_side = 2*sqrt(2)*r_outer_m;
    A_uc = L_side^2;
    flux_tot = flux_tot_per_m2*A_uc;

    eff_tot = Gtot/flux_tot;
    eff_vec_marginal = Gtot_vec/flux_tot;
    eff_vec = cumsum(Gtot_vec)/flux_tot;

    plot(1:num_collectors, eff_vec, colors{zind})
end

%%
xlim([1 num_collectors])
xlabel('Number of collection tiers')
ylabel('Collection Efficiency')
%ylim([0.05 0.30])
grid on
fixfigs(1,3,14,12)


%% efficiency vs radius for one and two tiers
num_collectors = 2;
r_col_vec = linspace(1e-6,20e-6,5);
num_radii = length(r_col_vec);


zo_m_vec = (0:4) * 1e-6;
figure(1)
clf
hold on
colors = {'k', 'b', 'g', 'y', 'r' };

Gtot_mat = zeros(num_collectors, num_radii);
for nind = 1:num_collectors
    
    for rind = 1:num_radii
        r_outer_m = r_col_vec(rind);
        depth_m = 0.5*r_outer_m  + (nind-1)*r_outer_m*sqrt(2);
        
        [G_shell, G_inner, G_outer, flux_tot, eff_col_shell, flux_tot_per_m2, i_act, J_act] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, depth_m);
        
        Gtot_mat(nind,rind) = 2*G_shell; % fcc lattice -- two shells per layer

    end
end

L_side = 2*sqrt(2)*r_outer_m;
A_uc = L_side^2;
flux_tot = flux_tot_per_m2*A_uc;

eff_mat_marginal = Gtot_mat/flux_tot;
eff_mat = cumsum(Gtot_mat)/flux_tot;

figure(2)
clf
hold on
plot(r_col_vec*1e6, eff_mat(1,:), 'b')
plot(r_col_vec*1e6, eff_mat(2,:), 'r')
xlabel('Collection Radius (um)')
ylabel('Collection Efficiency')
grid on
fixfigs(2,3,14,12)
