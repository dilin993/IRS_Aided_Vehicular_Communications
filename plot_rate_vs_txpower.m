%% IRS aidded V2X Uplink Scenario

clear all; 
close all; clc;

gamma_th_dB = 60; % SNR threshold for outage
sigma_n_sqr_dBm = -60; % Noise Power
p_sig_dBm = 0:0.5:30; % Signal power
n_trials = 10; % no. of trials
n_angle_realizations = 20; % angle realizations
fc = 24.2; % 24.2 GHz
error_th = 1e-3;
path_loss_enabled = 1;
K_r = 1; % IRS subgroup row size
K_c = 1; % IRS subgroup column size
sigma_sqr = 1; % Variance of Rician distribution
L=8; % number of discrete levels
early_stop = 1;
iteration_count = 100;
random_initilization = 0;

sim1 = simulation(sigma_n_sqr_dBm, gamma_th_dB,0, n_trials,...
    fc, error_th, n_angle_realizations, path_loss_enabled, K_r,...
    K_c, sigma_sqr,iteration_count, L, early_stop,random_initilization);
K_r = 2;
K_c = 2;
sim2 = simulation(sigma_n_sqr_dBm, gamma_th_dB,0, n_trials,...
    fc, error_th, n_angle_realizations, path_loss_enabled, K_r,...
    K_c, sigma_sqr,iteration_count, L, early_stop,random_initilization);


% define IRS
id_irs = 1;
pos_irs = [0;0;1];
Nx_irs = 1;
Ny_irs = 16;
Nz_irs = 16;
delta_irs = 0.5;
gain_irs = 10;
irs_256 = node(id_irs, Nx_irs, Ny_irs, Nz_irs, delta_irs, pos_irs, gain_irs);
Ny_irs = 8;
Nz_irs = 8;
irs_64 = node(id_irs, Nx_irs, Ny_irs, Nz_irs, delta_irs, pos_irs, gain_irs);

%define BS
id_bs = 2;
pos_bs = [20;-10;2];
Nx_bs = 4;
Ny_bs = 1;
Nz_bs = 2;
delta_bs = 0.5;
gain_bs = 5;
bs = node(id_bs, Nx_bs, Ny_bs, Nz_bs, delta_bs, pos_bs, gain_bs);

% define car
id_car = 2;
pos_car = [1.5;0;1];
Nx_car = 1;
Ny_car = 1;
Nz_car = 1;
delta_car = 0.5;
gain_car = 5;
car = node(id_car, Nx_car, Ny_car, Nz_car, delta_car, pos_car, gain_car);

num_points = size(p_sig_dBm, 2);

R_irs_256_full = zeros(num_points,1);
R_irs_64_full = zeros(num_points,1);
R_irs_256 = zeros(num_points,1);
R_irs_64 = zeros(num_points,1);
R_irs_256_full_sub = zeros(num_points,1);
R_irs_64_full_sub = zeros(num_points,1);
R_irs_256_sub = zeros(num_points,1);
R_irs_64_sub = zeros(num_points,1);
R_nirs = zeros(num_points,1);


parfor j=1:sim1.angle_realizations
    fprintf('Iterating angle realization %d\n', j);
    for n=1:sim1.n_trials
        [R_irs_256_1_full, R_irs_256_1, R_nirs_1] = simulate(...
            bs,irs_256,car, sim1,p_sig_dBm);
        [R_irs_64_1_full,R_irs_64_1, ~] = simulate(bs, irs_64,...
            car, sim1, p_sig_dBm);
        
        [R_irs_256_1_full_sub, R_irs_256_1_sub, ~] = ...
            simulate(bs,irs_256,car, sim2,p_sig_dBm);
        [R_irs_64_1_full_sub,R_irs_64_1_sub, ~] = simulate(bs,...
            irs_64,car, sim2, p_sig_dBm);
        
        
        R_irs_256_full = R_irs_256_full + R_irs_256_1_full;
        R_irs_64_full = R_irs_64_full + R_irs_64_1_full;
        R_irs_256 = R_irs_256 + R_irs_256_1;
        R_irs_64 = R_irs_64 + R_irs_64_1;
        R_irs_256_full_sub = R_irs_256_full_sub + R_irs_256_1_full_sub;
        R_irs_64_full_sub = R_irs_64_full_sub + R_irs_64_1_full_sub;
        R_irs_256_sub = R_irs_256_sub + R_irs_256_1_sub;
        R_irs_64_sub = R_irs_64_sub + R_irs_64_1_sub;
        R_nirs = R_nirs + R_nirs_1;
    end
    irs_256.change_angles();
    irs_64.change_angles();
    bs.change_angles();
    car.change_angles();
end

%%
R_irs_256_full = R_irs_256_full/(sim1.n_trials*sim1.angle_realizations);
R_irs_64_full = R_irs_64_full/(sim1.n_trials*sim1.angle_realizations);
R_irs_256 = R_irs_256/(sim1.n_trials*sim1.angle_realizations);
R_irs_64 = R_irs_64/(sim1.n_trials*sim1.angle_realizations);
R_irs_256_full_sub = R_irs_256_full_sub/(sim1.n_trials*sim1.angle_realizations);
R_irs_64_full_sub = R_irs_64_full_sub/(sim1.n_trials*sim1.angle_realizations);
R_irs_256_sub = R_irs_256_sub/(sim1.n_trials*sim1.angle_realizations);
R_irs_64_sub = R_irs_64_sub/(sim1.n_trials*sim1.angle_realizations);
R_nirs = R_nirs/(sim1.n_trials*sim1.angle_realizations);
save('data_txpower3');

%%
fig1 = figure;
hold on;
plot(p_sig_dBm, R_irs_256_full,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64_full,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_256,'--','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64,'--','LineWidth',1.5);
plot(p_sig_dBm, R_irs_256_full_sub,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64_full_sub,'-','LineWidth',1.5);
% plot(p_sig_dBm, R_irs_256_sub,'--','LineWidth',1.5);
% plot(p_sig_dBm, R_irs_64_sub,'--','LineWidth',1.5);
plot(p_sig_dBm, R_nirs,'--','LineWidth',1.5);
hold off;
grid on;
xlabel('Tx Power(dBm)');
ylabel('Achievable Rate (bits/s/Hz)');
% legend('with IRS(256) full CSI','with IRS(64) full CSI', ...
%     'with IRS(256) LoS','with IRS(64) LoS',...
%     'with IRS(256) full CSI subgroups','with IRS(64) full CSI subgroups', ...
%     'with IRS(256) LoS subgroups','with IRS(64) LoS subgroups',...
%     'without IRS','Location','best');
legend('with IRS(256) full CSI','with IRS(64) full CSI', ...
       'with IRS(256) position based','with IRS(64) position based',...
       'with IRS(256) 2x2 grouping','with IRS(64) 2x2 grouping',...
       'without IRS','Location','best');
title('Achievable Rate with Transmit Power');
% date = datestr(now,'YYYY.mm.dd.HH.MM');
print(gcf,'achievable_rate_vs_txpower_all.png','-dpng','-r400');

%%
function [R_irs_full,R_irs,R_nirs] = ...
    simulate(bs, irs, car,sim, p_sig_dBm)
    R_irs_full = zeros(size(p_sig_dBm,2),1);
    R_irs = zeros(size(p_sig_dBm,2),1);
    R_nirs = zeros(size(p_sig_dBm,2),1);
    d_d = dist_3d(bs.pos, car.pos);
    d_r = dist_3d(bs.pos, irs.pos);
    d_v = dist_3d(irs.pos, car.pos);
    h_d = generate_MIMO_channel(car, bs, d_d,...
        sim.fc, 0, 1e15, 1, sim.path_loss_enabled);
    H_r = generate_MIMO_channel(irs, bs, d_r,...
        sim.fc, 0, 2, 2, sim.path_loss_enabled);
    h_v = generate_MIMO_channel(car, irs, d_v,...
        sim.fc, 0, 1, 1, sim.path_loss_enabled);
    h_d_1 = generate_MIMO_channel(car, bs, d_d,...
        sim.fc, 0, 0, 1, sim.path_loss_enabled);
    H_r_1 = generate_MIMO_channel(irs, bs, d_r,...
        sim.fc, 0, 0, 2, sim.path_loss_enabled);
    h_v_1 = generate_MIMO_channel(car, irs, d_v,...
        sim.fc, 0, 0, 1, sim.path_loss_enabled);
    for i=1:size(p_sig_dBm,2)
        p_sig = db2pow(p_sig_dBm(i)-30);
        % passive beamforming with full CSI
        [~,R_irs_full(i)] = optimal_phase_shift3(h_d, H_r, h_v,...
            sim.error_th,sim.sigma_n_sqr,p_sig, sim.K_r, sim.K_c,...
            irs.Ny, irs.Nz,sim.iteration_count, sim.L, sim.early_stop,...
            sim.random_initilization);
        R_nirs(i)=achievable_rate(h_d, sim.sigma_n_sqr,p_sig);
        
         % position based passive beamforming
        [v,~] = optimal_phase_shift3(h_d_1, H_r_1, h_v_1,...
            sim.error_th,sim.sigma_n_sqr, p_sig, sim.K_r, sim.K_c,...
            irs.Ny, irs.Nz,sim.iteration_count, sim.L, sim.early_stop,...
            sim.random_initilization);
        R_irs(i)=achievable_rate(h_d + H_r*diag(v)*h_v,...
            sim.sigma_n_sqr, p_sig);
    end
end










