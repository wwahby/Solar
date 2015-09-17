function Jsc_z = calc_Jsc_z(z, m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot)
    lambda_vec = linspace(lambda_min, lambda_max, Npts_lambda);

    Jsc_lambda_z_vec = zeros(1, Npts_lambda);
    for lind = 1:Npts_lambda
        lambda = lambda_vec(lind);
        flux_per_m2 = vsc.get_photon_flux_per_m2( lambda );
        N = flux_per_m2/100^2; % N is per cm^2
        alpha = vsc.get_absorption_coefficient(lambda);
        
        Jsc_lambda_z_vec(lind) = vsc.calc_Jsc_lambda_z( z, m_max, W, d, Dn, tau_n, Sn_top, Sn_bot, alpha, N);
    end
    
    Jsc_z = trapz(lambda_vec, Jsc_lambda_z_vec);
    
end