lambda_nm_vec = linspace(1,1100, 1e2);


r_m = 1e-6;
zo_m = 20e-6;

nvol_vec = zeros(1,lind);
for lind = 1:length(lambda_nm_vec)
    lambda_nm = lambda_nm_vec(lind);
    nvol = calc_nvol(lambda_nm, r_m, zo_m);
    nvol_vec(lind) = nvol;
end

nvol_tot = trapz(lambda_nm_vec, nvol_vec);

figure(1)
clf
semilogy(lambda_nm_vec, nvol_vec)
ylim([1e-10 1e10])
fixfigs(1,3,14,12)

figure(2)
clf
plot(lambda_nm_vec, nvol_vec)
fixfigs(2,3,14,12)