function n = calc_n(Npts_x, Npts_z, W, d, N, Dn, tau_n, alpha, Sn_top, Sn_bot, m_max)
    % Npts_x - (-) Number of points to use for horizontal dimension
    % Npts_z - (-) Number of points to use for vertical dimension
    % zvec - (cm) vector of z positions to calculate z for
    % W - (cm) width of central absorption region
    % d - (cm) depth of central absorption region (should be positive)
    % N - (#/cm^2) areal photon density at wavelength of interest
    % Dn - (cm^2/s) Diffusion coefficient of electrons in the p-doped central region
    % tau_n - (s) Electron lifetime in the p-doped central region
    % alpha - (-) absorption coefficient of the wavelength under consideration
    % Sn_top - (cm/s) recombination velocity at the top surface
    % Sn_bot - (cm/s) recombination velocity at the bottom surface
    % m_max - (-) maximum number of terms of the fourier series to include
    
    % position vectors
    xvec = linspace(-W/2, W/2, Npts_x);
    zvec = linspace(0, d, Npts_z);
    % Diffusion Length
    Ln = sqrt(Dn * tau_n);
    
    % Precalculate coefficients
    [A_m_vec, B_m_vec, beta_m_vec] = vsc.calc_vert_cell_coefficients(m_max,  W, d, Dn, Ln, Sn_top, Sn_bot, alpha );
    
    % Calculate n everywhere
    n = zeros( length(xvec), length(zvec) );
    for xind = 1:length(xvec)
        x = xvec(xind);
        for zind = 1:length(zvec)
            z = zvec(zind);
            
            nxz = 0;
            for m = 1:m_max
                nxz = nxz + vsc.calc_n_m(m, x, z, W, N, Ln, alpha, beta_m_vec(m), A_m_vec(m), B_m_vec(m));
            end
            n(xind, zind) = nxz;
        end
    end
    
end