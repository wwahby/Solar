function alpha_cm = get_absorption_coefficient(lambda_nm)
    [lambda_nm_vec, absorption_vec_cm, n, k] = vsc.get_si_optical_data();
    
    if (lambda_nm < min(lambda_nm_vec))
        alpha_cm = max(absorption_vec_cm);
    elseif (lambda_nm > max(lambda_nm_vec))
        alpha_cm = min(absorption_vec_cm);
    else
        alpha_cm = interp1(lambda_nm_vec, absorption_vec_cm, lambda_nm, 'pchip', 'extrap');
    end
end