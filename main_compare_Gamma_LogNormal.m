clear all
close all
%
snrdB = -10:1:15; % average transmit SNR

% figure(1);
% load('data_Gamma_30_25_20.mat')
% semilogy(snrdB,P_out_CCT_Gamma_ana,'r-','linewidth',2); hold on;
% load('data_LogNormal_30_25_20.mat')
% semilogy(snrdB,P_out_CCT_LN_ana,'b--','linewidth',2); hold on;
% axis([-10 5 10^(-6) 10^(0)]);
% legend('using Gamma',...
%     'using LogNormal',...
%     'Location','sw');
% 
% load('data_Gamma_30_25_20.mat')
% zoomPlot(snrdB,P_out_CCT_Gamma_sim,[-5 -4],[0.6 0.2 0.2 0.7],[0 0]);
% hold on
% load('data_LogNormal_30_25_20.mat')
% semilogy(snrdB,P_out_CCT_LN_ana,'b--','linewidth',2); hold on;
% legend hide


%-----Plot EC
fig2 = figure(2);

% for CCT
load('data_Gamma_30_25_20.mat')
% plot(snrdB,EC_CCT_sim,'r^'); hold on;
plot(snrdB,EC_CCT_ana,'r-','linewidth',2); hold on;

% for OST
load('data_LogNormal_30_25_20.mat')
% plot(snrdB,EC_CCT_sim,'b^'); hold on;
plot(snrdB,EC_CCT_ana,'b--','linewidth',2); hold on;

xlabel('Transmit SNR (dB)');
ylabel('Ergodic Capacity (bits/sec/Hz');
legend('using Gamma distribution',...
    'using LogNormal distribution',...
    'Location','se');

load('data_Gamma_30_25_20.mat')
zoomPlot(snrdB,EC_CCT_ana,[11 12],[0.2 0.55 0.33 0.3],[0 0]);
% plot(snrdB,EC_CCT_ana,'r-','linewidth',2); hold on;
hold on
load('data_LogNormal_30_25_20.mat')
plot(snrdB,EC_CCT_ana,'b--','linewidth',2); hold on;
legend hide
