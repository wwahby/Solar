function [eta_col_vec, eta_abs_vec, flux_vec_m2, alpha_vec_cm, flux_wehrli_vec_m2] = calc_collection_efficiency(lambda_vec_nm, m_max, W_cm, d_cm, Dn_cm2, Ln_cm, Sn_cm)

    eta_col_vec = zeros(1, length(lambda_vec_nm));
    eta_abs_vec = zeros(1, length(lambda_vec_nm));
    flux_vec_m2 = zeros(1, length(lambda_vec_nm));
    flux_wehrli_vec_m2 = zeros(1, length(lambda_vec_nm));
    alpha_vec_cm = zeros(1, length(lambda_vec_nm));
    for lind = 1:length(lambda_vec_nm)
        lambda_nm = lambda_vec_nm(lind);
        
        [flux_per_m2, flux_wehrli_per_m2] = vsc.get_photon_flux_per_m2( lambda_nm );
        N_cm2 = flux_per_m2/100^2; % N needs to be in cm^-2
        N_cm2_wehrli = flux_wehrli_per_m2/1e4;
        
        alpha_cm = vsc.get_absorption_coefficient(lambda_nm);
        alpha_vec_cm(lind) = alpha_cm;
        flux_vec_m2(lind) = flux_per_m2;
        flux_wehrli_vec_m2(lind) = flux_wehrli_per_m2;
        
        [eta_col_lambda, eta_abs_lambda] = vsc.calc_collection_efficiency_lambda(m_max, W_cm, d_cm, Dn_cm2, Ln_cm, Sn_cm, alpha_cm, N_cm2);
        eta_col_vec(lind) = eta_col_lambda;
        eta_abs_vec(lind) = eta_abs_lambda;
    end

end