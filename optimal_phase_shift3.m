function [v, R] = optimal_phase_shift3(h_d, H_r, h_v, error_th,...
    sigma_n_sqr, p_sig, K_r, K_c, N_r, N_c, iteration_count, L,...
    early_stop,random_initilization)

if K_r==1 && K_c==1
    [v,R] = optimal_phase_shift(h_d, H_r, h_v, error_th, sigma_n_sqr,...
    p_sig, iteration_count, L, early_stop,random_initilization);
    return;
end

H_r_exp = zeros(N_r, N_c, size(H_r,1));
for i=1:size(H_r,1)
    H_r_exp(:,:,i) = reshape(H_r(i,:),N_r,N_c);
end
h_v_exp = reshape(h_v,N_r,N_c);
H_r1 = zeros(size(H_r,1),size(H_r,2)/(K_r*K_c));
h_v1 = reshape(h_v_exp(1:K_r:end,1:K_c:end),N_r*N_c/(K_r*K_c),1);
for i=1:size(H_r,1)
    H_r1(i,:) = reshape(H_r_exp(1:K_r:end,1:K_c:end,i),N_r*N_c/(K_r*K_c),1);
end

[v1,~,~] = optimal_phase_shift(h_d, H_r1, h_v1, error_th, sigma_n_sqr,...
    p_sig, iteration_count, L, early_stop,random_initilization);
v1_exp = reshape(v1,N_r/K_r,N_c/K_c);
v = reshape(imresize(v1_exp,[N_r,N_c]),N_r*N_c,1);
v = exp(1i*angle(v));
R = achievable_rate(h_d + H_r*diag(v)*h_v, sigma_n_sqr,p_sig);
end