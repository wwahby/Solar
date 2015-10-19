function [eta_col_lambda, eta_abs_lambda] = calc_collection_efficiency_lambda(m_max, W, d_cm, Dn, Ln, Sn, alpha_cm, N)
    q = 1.602e-19; % (C) electron charge
    
    eta_abs_lambda = 1 - exp(-alpha_cm*d_cm); % absorption efficiency

    Isc_lambda = vsc.calc_Isc_lambda(m_max, W, d_cm, Dn, Ln, Sn, alpha_cm, N);
    eta_col_lambda = Isc_lambda/(W*q*N*eta_abs_lambda);
end