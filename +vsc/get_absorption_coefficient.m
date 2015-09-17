function alpha = get_absorption_coefficient(lambda)
    [lambda_nm_vec, absorption_vec, n, k] = vsc.get_si_optical_data();
    alpha = interp1(lambda_nm_vec, absorption_vec, lambda, 'pchip', 'extrap');
end