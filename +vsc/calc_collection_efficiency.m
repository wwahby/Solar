function [eta_col_vec, eta_abs_vec, flux_vec, alpha_vec] = calc_collection_efficiency(lambda_vec, m_max, W, d, Dn, Ln, Sn)

    eta_col_vec = zeros(1, length(lambda_vec));
    eta_abs_vec = zeros(1, length(lambda_vec));
    flux_vec = zeros(1, length(lambda_vec));
    alpha_vec = zeros(1, length(lambda_vec));
    for lind = 1:length(lambda_vec)
        lambda = lambda_vec(lind);
        
        flux_per_m2 = vsc.get_photon_flux_per_m2( lambda );
        N = flux_per_m2/100^2; % N needs to be in cm^-2
        alpha = vsc.get_absorption_coefficient(lambda);
        alpha_vec(lind) = alpha;
        flux_vec(lind) = flux_per_m2;
        
        [eta_col_lambda, eta_abs_lambda] = vsc.calc_collection_efficiency_lambda(m_max, W, d, Dn, Ln, Sn, alpha, N);
        eta_col_vec(lind) = eta_col_lambda;
        eta_abs_vec(lind) = eta_abs_lambda;
    end

end