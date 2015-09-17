function Isc_lambda_m = calc_Isc_lambda_m(W, d, N, alpha, beta_m, Dn, Sn)
    % Assumes Sn_bot = Sn_top = Sn
    q = 1.602e-19; % (C) electron charge
    Isc_lambda_m = 8*N*q/W * (1-exp(-alpha*d))/(beta_m^2 - alpha^2) ...
        * ( 1 - alpha/beta_m * (Sn * coth(1/2*alpha*d) + alpha*Dn)/(Sn*coth(1/2*beta_m*d) + beta_m*Dn) );
end