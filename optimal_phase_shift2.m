function [v, R] = optimal_phase_shift2(h_d, H_r, h_v, error_th,...
    sigma_n_sqr, p_sig, K_r, K_c, N_r, N_c, iteration_count, L, early_stop)
H_r1 = zeros(size(H_r,1),size(H_r,2)/(K_r*K_c));
h_v1 = zeros(size(h_v,1)/(K_r*K_c),1);
n=1;
for i=1:K_r:N_r
    for j=1:K_c:N_c
        H_r1(:,n) = H_r(:,(i-1)*N_c + j);
        h_v1(n) = h_v((i-1)*N_c + j);
        n = n + 1;
    end
end
[v1,~,~] = optimal_phase_shift(h_d, H_r1, h_v1, error_th, sigma_n_sqr,...
    p_sig, iteration_count, L, early_stop);
v = zeros(N_r*N_c,1);
n=1;
for i=1:K_r:N_r
    for j=1:K_c:N_c
        for l=1:K_r
            for m=1:K_c
                v((i-1+l-1)*N_c + j+m-1) = v1(n);
            end
        end
        n = n + 1;
    end
end
R = achievable_rate(h_d + H_r*diag(v)*h_v, sigma_n_sqr,p_sig);
end