clear all; 
close all; clc;

load('data_txpower2');

fig1 = figure;
hold on;
plot(p_sig_dBm, R_irs_256_full,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64_full,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_256,'--','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64,'--','LineWidth',1.5);
plot(p_sig_dBm, R_irs_256_full_sub,'--','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64_full_sub,'--','LineWidth',1.5);
plot(p_sig_dBm, R_nirs,'-','LineWidth',1.5);
hold off;
grid on;
xlabel('Tx Power (dBm)');
ylabel('Achievable Rate (bits/s/Hz)');
legend('with IRS(256) full CSI','with IRS(64) full CSI', ...
       'with IRS(256) position based','with IRS(64) position based',...
       'with IRS(256) 2x2 grouping','with IRS(64) 2x2 grouping',...
       'without IRS','Location','best');
title('Full CSI vs. Position based Beamforming');
print(gcf,'achievable_rate_vs_txpower_pos.png','-dpng','-r400');


fig2 = figure;
hold on;
plot(p_sig_dBm, R_irs_256_full,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64_full,'-','LineWidth',1.5);
plot(p_sig_dBm, R_irs_256_full_sub,'--','LineWidth',1.5);
plot(p_sig_dBm, R_irs_64_full_sub,'--','LineWidth',1.5);
plot(p_sig_dBm, R_nirs,'-','LineWidth',1.5);
hold off;
grid on;
xlabel('Tx Power (dBm)');
ylabel('Achievable Rate (bits/s/Hz)');
legend('with IRS(256) full CSI','with IRS(64) full CSI', ...
   'with IRS(256) 2x2 grouping','with IRS(64) 2x2 grouping', ...
       'without IRS','Location','best');
title('Full CSI vs. Grouping based Beamforming');
print(gcf,'achievable_rate_vs_txpower_sub.png','-dpng','-r400');


