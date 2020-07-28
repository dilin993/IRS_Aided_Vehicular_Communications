function [v, R, history] = optimal_phase_shift(h_d, H_r, h_v, error_th,...
    sigma_n_sqr, p_sig, iteration_count, L, early_stop,random_initilize)
Phi = H_r * diag(h_v);
A = Phi'*Phi;
b = Phi'*h_d;
discrete_levels = linspace(0, 2*pi, L);
N = size(H_r,2);
if random_initilize==0
    I = ones(N,1);
else
    I = randi([1 L],N,1);
end
for n=1:N
    theta(n) = discrete_levels(I(n));
end
v = exp(1i * theta);
obj_value = achievable_rate(h_d + H_r*diag(v)*h_v, sigma_n_sqr,p_sig);
history = zeros(iteration_count+1,1);
history(1,1) = obj_value;
for iter = 1:iteration_count
    for n=1:N
        kappa_n = b(n);
        for l=1:N
            if l~=n
                kappa_n = kappa_n + A(n,l)*v(l);
            end
        end
        min_diff = 4*pi;
        for ang=discrete_levels
            diff = abs(angle(exp(1i*ang))-angle(kappa_n));
            if diff < min_diff
                min_diff = diff;
                theta(n) = ang;
            end
        end
    end
    v = exp(1i * theta);
    prev_obj_value = obj_value;
    obj_value=achievable_rate(h_d + H_r*diag(v)*h_v, sigma_n_sqr,p_sig);
    history(iter + 1,1) = obj_value;
    if early_stop==1 && abs(obj_value - prev_obj_value) < error_th
        break;
    end
end
R = obj_value;
v = exp(1i * theta);
end