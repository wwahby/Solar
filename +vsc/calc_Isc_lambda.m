function Isc_lambda = calc_Isc_lambda(m_max, W, d, Dn, Ln, Sn, alpha, N)
    % Assumes Sn_bot = Sn_top = Sn

    % Precalculate coefficients
    [A_m_vec, B_m_vec, beta_m_vec] = vsc.calc_vert_cell_coefficients(m_max,  W, d, Dn, Ln, Sn, Sn, alpha );
    
    Isc_lambda = 0;
    for m = 1:m_max
        Isc_lambda = Isc_lambda + vsc.calc_Isc_lambda_m(W, d, N, alpha, beta_m_vec(m), Dn, Sn);
    end
    
end