function [eta_col_lambda, eta_abs_lambda] = calc_collection_efficiency_lambda(m_max, W, d, Dn, Ln, Sn, alpha, N)
    q = 1.602e-19; % (C) electron charge
    
    eta_abs_lambda = 1 - exp(-alpha*d); % absorption efficiency

    Isc_lambda = vsc.calc_Isc_lambda(m_max, W, d, Dn, Ln, Sn, alpha, N);
    eta_col_lambda = Isc_lambda/(W*q*N*eta_abs_lambda);
end