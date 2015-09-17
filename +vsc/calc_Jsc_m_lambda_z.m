function Jsc_m_lambda_z = calc_Jsc_m_lambda_z(z, W, alpha, beta_m, A_m, B_m, N)
    q = 1.602e-19; % (C) electron charge
    Jsc_m_lambda_z = 8*N*alpha*q/W/(beta_m^2 - alpha^2) * (exp(-alpha*z) + A_m*exp(beta_m*z) + B_m*exp(-beta_m*z) );
end