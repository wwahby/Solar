function Jsc_lambda_z = calc_Jsc_lambda_z( z, m_max, W, d, Dn, tau_n, Sn_top, Sn_bot, alpha, N)
    % Calculates short circuit current density at a position z and a wavelength
    % lambda (implicit in N)
    % z - (cm) vertical coordinate under consideration
    % W - (cm) width of central absorption region
    % d - (cm) depth of central absorption region
    % N - (#/cm^2) areal photon density at wavelength of interest
    % Dn - (cm^2/s) Diffusion coefficient of electrons in the p-doped central region
    % tau_n - (s) Electron lifetime in the p-doped central region
    % alpha - (-) absorption coefficient of the wavelength under consideration
    % Sn_top - (cm/s) recombination velocity at the top surface
    % Sn_bot - (cm/s) recombination velocity at the bottom surface
    % m_max - (-) maximum number of terms of the fourier series to include

    % Diffusion Length
    Ln = sqrt(Dn * tau_n);
    
    % Precalculate coefficients
    [A_m_vec, B_m_vec, beta_m_vec] = vsc.calc_vert_cell_coefficients(m_max,  W, d, Dn, Ln, Sn_top, Sn_bot, alpha );
    
    Jsc = 0;
    for m = 1:m_max
        Jsc_new_m = vsc.calc_Jsc_m_lambda_z(z, W, alpha, beta_m_vec(m), A_m_vec(m), B_m_vec(m), N);
        Jsc_new_m( isnan(Jsc_new_m) ) = 0; % clear NaNs
        Jsc = Jsc + Jsc_new_m; %vsc.calc_Jsc_m_lambda_z(z, W, alpha, beta_m_vec(m), A_m_vec(m), B_m_vec(m), N);
    end
    Jsc_lambda_z = Jsc;
end