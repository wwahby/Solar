lambda_nm_vec = linspace(1,1100, 1e2);

r_dep_m = 2e-6;
L_col_m = 10e-6;
r_col_m = r_dep_m + L_col_m;
zo_m = 12e-6;

r_inner_m = r_dep_m;
r_outer_m = r_col_m;
zo_top_m = zo_m;
zo_bot_m = zo_top_m + sqrt(3)*r_outer_m;
[G_shell_top, G_inner_top, G_outer_top, flux_tot, eff_col_shell_top, i_act_top, J_act_top] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_top_m);
[G_shell_bot, G_inner_bot, G_outer_bot, flux_tot, eff_col_shell_bot, i_act_bot, J_act_bot] = calc_generation_rate_in_spherical_shell(lambda_nm_vec, r_inner_m, r_outer_m, zo_bot_m);

G_tot = G_shell_top + G_shell_bot
i_act_tot = i_act_top + i_act_bot
J_act_tot = J_act_top + J_act_bot

eff_tot = G_tot/flux_tot