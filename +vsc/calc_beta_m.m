function beta_m = calc_beta_m( m, Ln, W)

    beta_m = sqrt(1/Ln^2 + (2*m-1)^2 * pi^2/W^2);
end