function [eta_col_vec, eta_abs_vec, flux_vec_m2, alpha_vec_cm] = calc_collection_efficiency(lambda_vec_nm, m_max, W, d, Dn, Ln, Sn)

    eta_col_vec = zeros(1, length(lambda_vec_nm));
    eta_abs_vec = zeros(1, length(lambda_vec_nm));
    flux_vec_m2 = zeros(1, length(lambda_vec_nm));
    alpha_vec_cm = zeros(1, length(lambda_vec_nm));
    for lind = 1:length(lambda_vec_nm)
        lambda_nm = lambda_vec_nm(lind);
        
        flux_per_m2 = vsc.get_photon_flux_per_m2( lambda_nm );
        N = flux_per_m2/100^2; % N needs to be in cm^-2
        alpha_cm = vsc.get_absorption_coefficient(lambda_nm);
        alpha_vec_cm(lind) = alpha_cm;
        flux_vec_m2(lind) = flux_per_m2;
        
        [eta_col_lambda, eta_abs_lambda] = vsc.calc_collection_efficiency_lambda(m_max, W, d, Dn, Ln, Sn, alpha_cm, N);
        eta_col_vec(lind) = eta_col_lambda;
        eta_abs_vec(lind) = eta_abs_lambda;
    end

end