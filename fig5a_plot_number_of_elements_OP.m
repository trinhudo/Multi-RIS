clear;
close all;

P_S_dB = -5:.3:25; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

%% ERA

load('data_Gamma_setL2.mat')
semilogy(P_S_dB,OP_ERA_Gamma_sim,'r-.', 'linewidth', 1.5); hold on;

load('data_Gamma_setL3.mat')
semilogy(P_S_dB,OP_ERA_Gamma_sim,'b-', 'linewidth', 1.5); hold on;

load('data_Gamma_setL4.mat')
semilogy(P_S_dB,OP_ERA_Gamma_sim,'g--', 'linewidth', 1.5); hold on;

%% ORA

load('data_Gamma_setL2.mat')
semilogy(P_S_dB,OP_ORA_Gamma_sim,'r-.', 'linewidth', 1.5); hold on;

load('data_Gamma_setL3.mat')
semilogy(P_S_dB,OP_ORA_Gamma_sim,'b-', 'linewidth', 1.5); hold on;

load('data_Gamma_setL4.mat')
semilogy(P_S_dB,OP_ORA_Gamma_sim,'g--', 'linewidth', 1.5); hold on;

%% SISO

load('data_Gamma_setL2.mat')
semilogy(P_S_dB,OP_SISO,'k-', 'linewidth', 1.5); hold on;

%% Plot

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Outage probability', 'Interpreter', 'Latex');
legend('Setting $\bf\mathrm{L}_2$', ...
    'Setting $\bf\mathrm{L}_3$', ...
    'Setting $\bf\mathrm{L}_4$', ...
    'Location','se',...
    'Interpreter', 'Latex');
axis([-5 25 10^(-5) 1e0]);
xticks([-5 0 5 10 15 20 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

