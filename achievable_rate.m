function R = achievable_rate(h, sigma_n_sqr, p_sig)
    R = log2(1 + p_sig * norm(h)^2 /sigma_n_sqr);
end