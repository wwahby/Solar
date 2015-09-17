function [A_m_vec, B_m_vec, beta_m_vec] = calc_vert_cell_coefficients(m_max,  W, d, Dn, Ln, Sn_top, Sn_bot, alpha )

% Precalculate coefficients
    A_m_vec = zeros(1, m_max);
    B_m_vec = zeros(1, m_max);
    beta_m_vec = zeros(1, m_max);
    for m = 1:m_max
        beta_m_vec(m) = vsc.calc_beta_m( m, Ln, W);
        A_m_vec(m) = vsc.calc_A_m( Sn_top, Sn_bot, beta_m_vec(m), alpha, Dn, d);
        B_m_vec(m) = vsc.calc_B_m( Sn_top, Sn_bot, beta_m_vec(m), alpha, Dn, d);
    end

end