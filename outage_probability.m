function p_out = outage_probability(gain, L, sigma_sqr, sigma_n_sqr,...
    gamma_th,p_sig)
    lambda1 = mean(gain)/sigma_sqr;
    a = sqrt(2*lambda1);
    b = sqrt(2*sigma_n_sqr*gamma_th./(sigma_sqr*p_sig'));
    p_out = 1 - marcumq(a,b,L);
end