function Voc_z = calc_Voc_z(z, m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot, T, Jo)
    k = 1.3806e-23; % (J/K) Boltzmann Constant
    q = 1.602e-19; % (C) Electron Charge

    %Jsc_z = calc_Jsc_z(z, m_max, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot);
    Jsc_z = vsc.calc_Jsc_z(z, m_max, lambda_min, lambda_max, Npts_lambda, W, d, Dn, tau_n, Sn_top, Sn_bot);
    
    Voc_z = k*T/q * log( Jsc_z/Jo + 1);
end
