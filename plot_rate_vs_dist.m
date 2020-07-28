%% IRS aidded V2X Uplink Scenario

clear all; 
close all; clc;

gamma_th_dB = 60; % SNR threshold for outage
sigma_n_sqr_dBm = -60; % Noise Power
p_sig_dBm = 30; % Signal power
n_trials = 200; % no. of trials
n_angle_realizations = 200; % angle realizations
fc = 24.2; % 24.2 GHz
error_th = 1e-3;
path_loss_enabled = 1;
K_r = 1; % IRS subgroup row size
K_c = 1; % IRS subgroup column size
sigma_sqr = 0.5; % Variance of Rician distribution
L=8; % number of discrete levels
early_stop = 1;
iteration_count = 100;
random_initilization = 0;

sim = simulation(sigma_n_sqr_dBm, gamma_th_dB,p_sig_dBm, n_trials,...
    fc, error_th, n_angle_realizations, path_loss_enabled, K_r,...
    K_c, sigma_sqr,iteration_count, L,early_stop,random_initilization);

% define IRS
pos_irs = [0;0;1];

%define BS
pos_bs = [20;-10;2];


num_points = 100;
car_points = cell(num_points, 1);



R_irs_full = zeros(num_points,1);
R_irs = zeros(num_points,1);
R_nirs = zeros(num_points,1);
P_o_irs_full = zeros(num_points,1);
P_o_irs = zeros(num_points,1);
P_o_nirs = zeros(num_points,1);
c_v = zeros(num_points,1);

for j=1:sim.angle_realizations
    fprintf('Iterating angle realization %d\n', j);
    i=1;
    irs = node(1, 1, 16, 16, 0.5, pos_irs,10);
    bs = node(2, 4, 1, 2, 0.5, pos_bs,5);
    for y=linspace(-50,50,num_points)
        pos_car = [1.5;y;1];
        car_points{i} = node(3, 1, 1, 1, 0.5, pos_car,0);
        c_v(i) = y;
        i = i + 1;
    end
    parfor n=1:sim.n_trials
        [R_irs_1_full, R_irs_1, R_nirs_1, P_o_irs_1_full,  P_o_irs_1,...
            P_o_nirs_1] = simulate(bs, irs, car_points, sim);
        R_irs_full = R_irs_full + R_irs_1_full;
        R_irs = R_irs + R_irs_1;
        R_nirs = R_nirs + R_nirs_1;
        P_o_irs_full = P_o_irs_full + P_o_irs_1_full;
        P_o_irs = P_o_irs + P_o_irs_1;
        P_o_nirs = P_o_nirs + P_o_nirs_1;
    end
end

%%
fig1 = figure;
R_irs_full = R_irs_full/(sim.n_trials*sim.angle_realizations);
R_irs = R_irs/(sim.n_trials*sim.angle_realizations);
R_nirs = R_nirs/(sim.n_trials*sim.angle_realizations);
P_o_irs_full = P_o_irs_full/(sim.n_trials*sim.angle_realizations);
P_o_irs = P_o_irs/(sim.n_trials*sim.angle_realizations);
P_o_nirs = P_o_nirs/(sim.n_trials*sim.angle_realizations);
n=1:num_points;
hold on;
plot(c_v, R_irs_full,'-','LineWidth',1.5);
% plot(c_v, R_irs,'-','LineWidth',1.5);
plot(c_v, R_nirs,'--','LineWidth',1.5);
grid on;
% ylim([0 0.2]);
xlabel('c_v (m)');
ylabel('Achievable Rate (bits/s/Hz)');
legend('with IRS(256) full CSI','without IRS');
title('Achievable Rate at Different Tx Positions');
% date = datestr(now,'YYYY.mm.dd.HH.MM');
print(gcf,'achievable_rate_vs_dist.png','-dpng','-r400');
% saveas(fig1, strcat('achievable_rate_vs_dist','.png'));


%%
function [R_irs_full,R_irs, R_nirs, P_o_irs_full, P_o_irs, P_o_nirs] = ...
simulate(bs, irs, car_points, sim)
    R_irs_full = zeros(size(car_points,1),1);
    R_irs = zeros(size(car_points,1),1);
    R_nirs = zeros(size(car_points,1),1);
    P_o_irs_full = zeros(size(car_points,1),1);
    P_o_irs = zeros(size(car_points,1),1);
    P_o_nirs = zeros(size(car_points,1),1);
    for i=1:size(car_points,1)
        d_d = dist_3d(bs.pos, car_points{i}.pos);
        d_r = dist_3d(bs.pos, irs.pos);
        d_v = dist_3d(irs.pos, car_points{i}.pos);
        h_d = generate_MIMO_channel(car_points{i}, bs, d_d,...
            sim.fc, 0, 1e15, 1, sim.path_loss_enabled);
        H_r = generate_MIMO_channel(irs, bs, d_r,...
            sim.fc, 0, 2, 2, sim.path_loss_enabled);
        h_v = generate_MIMO_channel(car_points{i}, irs, d_v,...
            sim.fc, 0, 1, 1, sim.path_loss_enabled);

        % passive beamforming with full CSI
        [~,R_irs_full(i)] = optimal_phase_shift3(h_d, H_r, h_v,...
        sim.error_th,sim.sigma_n_sqr,sim.p_sig, sim.K_r, sim.K_c,...
        irs.Ny, irs.Nz, sim.iteration_count, sim.L, sim.early_stop,...
        sim.random_initilization);
        R_nirs(i)=achievable_rate(h_d, sim.sigma_n_sqr,sim.p_sig);
    end
end










