function Isc = calc_Isc(m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, Ln, Sn)
    lambda_vec = linspace(lambda_min, lambda_max, Npts_lambda);
    
	Isc_lambda_vec = zeros(1, Npts_lambda);
    for lind = 1:Npts_lambda
        lambda = lambda_vec(lind);
        lambda_nm = lambda;
        alpha = vsc.get_absorption_coefficient(lambda);
        flux_per_m2 = vsc.get_photon_flux_per_m2( lambda );
        N = flux_per_m2/100^2; % N is per cm^2
        
        Isc_lambda_vec(lind) = vsc.calc_Isc_lambda(m_max, W, d, Dn, Ln, Sn, alpha, N);
    end
    
    Isc = trapz(lambda_vec, Isc_lambda_vec);
end
