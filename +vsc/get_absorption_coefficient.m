function alpha_cm = get_absorption_coefficient(lambda_nm)
    [lambda_nm_vec, absorption_vec_cm, n, k] = vsc.get_si_optical_data();
    alpha_cm = interp1(lambda_nm_vec, absorption_vec_cm, lambda_nm, 'pchip', 'extrap');
end