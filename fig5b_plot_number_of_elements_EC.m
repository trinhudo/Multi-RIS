clear;
close all;

P_S_dB = -5:.3:25; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

%% ERA

load('data_Gamma_setL2.mat')
plot(P_S_dB,EC_ERA_Gamma_sim,'r-.', 'linewidth', 1.5); hold on;

load('data_Gamma_setL3.mat')
plot(P_S_dB,EC_ERA_Gamma_sim,'b-', 'linewidth', 1.5); hold on;

load('data_Gamma_setL4.mat')
plot(P_S_dB,EC_ERA_Gamma_sim,'g--', 'linewidth', 1.5); hold on;

%% ORA

load('data_Gamma_setL2.mat')
plot(P_S_dB,EC_ORA_Gamma_sim,'r-.', 'linewidth', 1.5); hold on;

load('data_Gamma_setL3.mat')
plot(P_S_dB,EC_ORA_Gamma_sim,'b-', 'linewidth', 1.5); hold on;

load('data_Gamma_setL4.mat')
plot(P_S_dB,EC_ORA_Gamma_sim,'g--', 'linewidth', 1.5); hold on;

%% SISO

load('data_Gamma_setL2.mat')
plot(P_S_dB,EC_SISO,'k-', 'linewidth', 1.5); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [b/s/Hz]', 'Interpreter', 'Latex');
legend('Setting $\bf\mathrm{L}_2$',...
    'Setting $\bf\mathrm{L}_3$',...
    'Setting $\bf\mathrm{L}_4$',...
    'Location','nw',...
    'Interpreter', 'Latex');
xticks([-5 0 5 10 15 20 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%--------------------------------------------------------------------------

load('data_Gamma_setL2.mat')
zoomPlot(P_S_dB,EC_ERA_Gamma_sim,[13, 16],[.5 .6 0.4 0.3],[0 0]); hold on
plot(P_S_dB,EC_ERA_Gamma_sim,'r-.', 'linewidth', 1.5); hold on;
load('data_Gamma_setL3.mat')
plot(P_S_dB,EC_ERA_Gamma_sim,'b-', 'linewidth', 1.5); hold on;
load('data_Gamma_setL4.mat')
plot(P_S_dB,EC_ERA_Gamma_sim,'g--', 'linewidth', 1.5); hold on;
plot(P_S_dB,EC_ORA_Gamma_sim,'g--', 'linewidth', 1.5); hold on;

legend hide