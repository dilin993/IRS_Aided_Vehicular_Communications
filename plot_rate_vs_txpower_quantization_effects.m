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

sim_8bit = simulation(sigma_n_sqr_dBm, gamma_th_dB,0, n_trials,...
    fc, error_th, n_angle_realizations, path_loss_enabled, K_r,...
    K_c, sigma_sqr,iteration_count, L, early_stop,random_initilization);
L=4;
sim_4bit = simulation(sigma_n_sqr_dBm, gamma_th_dB,0, n_trials,...
    fc, error_th, n_angle_realizations, path_loss_enabled, K_r,...
    K_c, sigma_sqr,iteration_count, L, early_stop,random_initilization);

L=2;
sim_2bit = simulation(sigma_n_sqr_dBm, gamma_th_dB,0, n_trials,...
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

R_irs_8bit = zeros(num_points,1);
R_irs_4bit = zeros(num_points,1);
R_irs_2bit = zeros(num_points,1);
R_nirs = zeros(num_points,1);

for j=1:sim_8bit.angle_realizations
    fprintf('Iterating angle realization %d\n', j);
    for n=1:sim_8bit.n_trials
        [R_irs_8bit_1, R_nirs_1] = simulate(...
            bs,irs_256,car, sim_8bit,p_sig_dBm);
        [R_irs_4bit_1, ~] = simulate(...
            bs,irs_256,car, sim_4bit,p_sig_dBm);
        [R_irs_2bit_1, ~] = simulate(...
            bs,irs_256,car, sim_2bit,p_sig_dBm);
        
        
        R_irs_8bit = R_irs_8bit + R_irs_8bit_1;
        R_irs_4bit = R_irs_4bit + R_irs_4bit_1;
        R_irs_2bit = R_irs_2bit + R_irs_2bit_1;
        R_nirs = R_nirs + R_nirs_1;
    end
    irs_256.change_angles();
    irs_64.change_angles();
    bs.change_angles();
    car.change_angles();
end

%%
R_irs_8bit = R_irs_8bit/(sim_8bit.n_trials*sim_8bit.angle_realizations);
R_irs_4bit = R_irs_4bit/(sim_8bit.n_trials*sim_8bit.angle_realizations);
R_irs_2bit = R_irs_2bit/(sim_8bit.n_trials*sim_8bit.angle_realizations);
R_nirs = R_nirs/(sim_8bit.n_trials*sim_8bit.angle_realizations);
save('data_txpower_bit');

%%
fig1 = figure;
hold on;
plot(p_sig_dBm, R_irs_8bit,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_4bit,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_2bit,'-','LineWidth',1.5);
plot(p_sig_dBm, R_nirs,'--','LineWidth',1.5);
hold off;
grid on;
xlabel('Tx Power (dBm)');
ylabel('Achievable Rate (bits/s/Hz)');
legend('3 bit phase shifts','2 bit phase shifts','1 bit phase shifts',...
       'without IRS','Location','best');
title('Effect of phase shift quantization');
print(gcf,'achievable_rate_vs_txpower_bit.png','-dpng','-r400');

%%
function [R_irs,R_nirs] = ...
    simulate(bs, irs, car,sim, p_sig_dBm)
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
    for i=1:size(p_sig_dBm,2)
        p_sig = db2pow(p_sig_dBm(i)-30);
        % passive beamforming with full CSI
        [~,R_irs(i)] = optimal_phase_shift3(h_d, H_r, h_v,...
            sim.error_th,sim.sigma_n_sqr,p_sig, sim.K_r, sim.K_c,...
            irs.Ny, irs.Nz,sim.iteration_count, sim.L, sim.early_stop,...
            sim.random_initilization);
        R_nirs(i)=achievable_rate(h_d, sim.sigma_n_sqr,p_sig);
    end
end










