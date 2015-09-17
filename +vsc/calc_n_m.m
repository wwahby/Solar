function n_m = calc_n_m(m, x, z, W, N, Ln, alpha, beta_m, A_m, B_m)
    n1 = 4*N/alpha*tau*(-1)^m/(pi*Ln^2*(2*m-1)*(alpha^2 - beta_m^2));
    n2 = exp(-alpha*z) + A_m*exp(beta_m*z) + B_m*exp(-beta_m*z);
    n3 = cos((2*m-1)*pi*x/W);
    
    n_m = n1*n2*n3;
end
