clear all; 
close all; clc;

p_sig_dBm = linspace(0, 30, 100);

d_256 = readtable(...
'C:\Users\dilin\Documents\Github\IRS_Learning\data_256.csv');
d_64 = readtable(...
'C:\Users\dilin\Documents\Github\IRS_Learning\data_64.csv');

fig1 = figure;
hold on;
plot(p_sig_dBm, d_256.capacityWithIrs,'-','LineWidth',1.5);
plot(p_sig_dBm, d_64.capacityWithIrs,'-','LineWidth',1.5);
plot(p_sig_dBm, d_256.capacityWithoutIrs,'--','LineWidth',1.5);
hold off;
grid on;
xlabel('Tx Power (dBm)');
ylabel('Achievable Rate (bits/s/Hz)');
legend('with IRS(256) full CSI','with IRS(64) full CSI','without IRS',...
'Location','best');


title('Achievable Rate in Ray Tracing Simulation');
% date = datestr(now,'YYYY.mm.dd.HH.MM');
% saveas(fig1, strcat('achievable_rate_vs_txpower','.png'));
print(gcf,'achievable_rate_vs_txpower_ray.png','-dpng','-r400');
