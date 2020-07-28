classdef simulation
    
    properties
        sigma_n_sqr;
        gamma_th;
        p_sig;
        n_trials;
        fc;
        error_th;
        angle_realizations;
        path_loss_enabled;
        K_r;
        K_c;
        sigma_sqr;
        iteration_count;
        L;
        early_stop;
        random_initilization;
    end
    
    methods
        function obj = simulation(sigma_n_sqr_dBm, gamma_th_dB,...
                p_sig_dBm, n_trials, fc, error_th, angle_realizations,...
                path_loss_enabled, K_r, K_c, sigma_sqr,iteration_count,...
                L, early_stop,random_initilization)
            obj.sigma_n_sqr = db2pow(sigma_n_sqr_dBm-30);
            obj.gamma_th = db2pow(gamma_th_dB);
            obj.p_sig = db2pow(p_sig_dBm-30);
            obj.n_trials = n_trials;
            obj.fc = fc;
            obj.error_th = error_th;
            obj.angle_realizations = angle_realizations;
            obj.path_loss_enabled = path_loss_enabled;
            obj.K_r = K_r;
            obj.K_c = K_c;
            obj.sigma_sqr = sigma_sqr;
            obj.iteration_count = iteration_count;
            obj.L = L;
            obj.early_stop = early_stop;
            obj.random_initilization = random_initilization;
        end
    end
end

