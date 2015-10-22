lambda_nm_vec = linspace(1,1100, 1e2);

r_dep_m = 1e-6;
L_col_m = 3e-6;
r_col_m = r_dep_m + L_col_m;
zo_m = 4e-6;

r_inner_m = r_dep_m;
r_outer_m = r_col_m;
zo_top_m = zo_m;
zo_bot_m = zo_top_m + sqrt(3)*r_outer_m;
[G_shell_top, G_inner_top, G_outer_top, flux_tot, eff_col_shell_top, i_act_top, J_act_top] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_top_m);
[G_shell_bot, G_inner_bot, G_outer_bot, flux_tot, eff_col_shell_bot, i_act_bot, J_act_bot] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_bot_m);

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
        depth_m = zo_m  + nind*r_outer_m*sqrt(3);
        [G_shell, G_inner, G_outer, flux_tot, eff_col_shell, i_act, J_act] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, depth_m);
        Gtot_vec(nind) = G_shell;
    end
    Gtot = sum(Gtot_vec);

    eff_tot = Gtot/flux_tot;
    eff_vec_marginal = Gtot_vec/flux_tot;
    eff_vec = cumsum(Gtot_vec)/flux_tot;

    plot(1:num_collectors, eff_vec, colors{zind})
end
xlim([1 num_collectors])
xlabel('Number of collectors')
ylabel('Collection Efficiency')
fixfigs(1,3,14,12)