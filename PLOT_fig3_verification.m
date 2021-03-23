
clear;
close all;

N_RIS = 5; %Number of RISs

P_S_dB = 0:30; %Transmit power of the source, dBm, e.g., 200mW = 23dBm
P_S_dB_EC = 0:3:30;
%% Plot OP

%
load('./data/data_Gamma_25x5.mat')

subplot(1,4,1);

semilogy(P_S_dB,OP_ERA_Gamma_sim,'r^'); hold on;
semilogy(P_S_dB,OP_ERA_Gamma_ana,'r-'); hold on;

semilogy(P_S_dB,OP_ORA_Gamma_sim,'bs'); hold on;
semilogy(P_S_dB,OP_ORA_Gamma_ana,'b-'); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Outage Probability', 'Interpreter', 'Latex');
legend('ERA (sim.)',...
    'ERA (ana. Eq.(xx))',...
    'ORA (sim.)',...
    'ORA (ana. Eq.(xx))',...
    'Location','ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

axis([0 30 10^(-5) 1e0]);

%--------------------------------------------------------------------------
subplot(1,4,2);

plot(P_S_dB,EC_ERA_Gamma_sim,'r^'); hold on;
plot(P_S_dB,EC_ERA_Gamma_ana,'r-'); hold on;

plot(P_S_dB,EC_ORA_Gamma_sim,'bs'); hold on;
plot(P_S_dB,EC_ORA_Gamma_ana,'b-'); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Ergodic Capacity (bit/sec/Hz)', 'Interpreter', 'Latex');
legend('ERA (sim.)',...
    'ERA (ana. Eq.(xx))',...
    'ORA (sim.)',...
    'ORA (ana. Eq.(xx))',...
    'Location','nw',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);
%--------------------------------------------------------------------------
load('./data/data_LogNormal_25x5.mat')

subplot(1,4,3);

semilogy(P_S_dB,OP_ERA_LN_sim,'r^'); hold on;
semilogy(P_S_dB,OP_ERA_LN_ana,'r-'); hold on;

semilogy(P_S_dB,OP_ORA_LN_sim,'bs'); hold on;
semilogy(P_S_dB,OP_ORA_LN_ana,'b-'); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Outage Probability', 'Interpreter', 'Latex');
legend('ERA (sim.)',...
    'ERA (ana. Eq.(xx))',...
    'ORA (sim.)',...
    'ORA (ana. Eq.(xx))',...
    'Location','ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%--------------------------------------------------------------------------
axis([0 30 10^(-5) 1e0]);

subplot(1,4,4);

plot(P_S_dB,EC_ERA_LN_sim,'r^'); hold on;
plot(P_S_dB,EC_ERA_LN_ana,'r-'); hold on;

plot(P_S_dB,EC_ORA_LN_sim,'bs'); hold on;
plot(P_S_dB,EC_ORA_LN_ana,'b-'); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Ergodic Capacity (bit/sec/Hz)', 'Interpreter', 'Latex');
legend('ERA (sim.)',...
    'ERA (ana. Eq.(xx))',...
    'ORA (sim.)',...
    'ORA (ana. Eq.(xx))',...
    'Location','nw',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%% Plot EC

% figure;
% % for ERA
% load('./data/data_Gamma_set1.mat')
% plot(P_S_dB, EC_SISO, 'md-'); hold on;
% 
% plot(P_S_dB,EC_ERA_Gamma_sim,'r^'); hold on;
% plot(P_S_dB,EC_ERA_Gamma_ana,'r-'); hold on;
% 
% load('./data/data_Gamma_set2.mat')
% plot(P_S_dB,EC_ERA_Gamma_sim,'bo'); hold on;
% plot(P_S_dB,EC_ERA_Gamma_ana,'b-'); hold on;
% 
% load('./data/data_Gamma_set3.mat')
% plot(P_S_dB,EC_ERA_Gamma_sim,'gs'); hold on;
% plot(P_S_dB,EC_ERA_Gamma_ana,'g-'); hold on;

% % for ORA
% load('./data/data_Gamma_set1.mat')
% plot(P_S_dB,OP_ORA_Gamma_sim,'r^'); hold on;
% plot(P_S_dB,OP_ORA_Gamma_ana,'r-'); hold on;
% 
% load('./data/data_Gamma_set2.mat')
% plot(P_S_dB,OP_ORA_Gamma_sim,'bo'); hold on;
% plot(P_S_dB,OP_ORA_Gamma_ana,'b-'); hold on;
% 
% load('./data/data_Gamma_set3.mat')
% plot(P_S_dB,OP_ORA_Gamma_sim,'gs'); hold on;
% plot(P_S_dB,OP_ORA_Gamma_ana,'g-'); hold on;

% xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
% ylabel('Ergodic Capacity', 'Interpreter', 'Latex');
% legend('SISO (sim.)',...
%     'Setting 1 (sim.)', 'Setting 1 (ana.)',...
%     'Setting 2 (sim.)', 'Setting 2 (ana.)',...
%     'Setting 3 (sim.)', 'Setting 3 (ana.)',...
%     'Location','ne',...
%     'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

% % for CCT
% load('./data/data_Gamma_25_25_25.mat')
% plot(P_S_dB,EC_ERA_sim,'r^'); hold on;
% plot(P_S_dB,EC_ERA_ana,'r-'); hold on;
% load('./data/data_Gamma_30_25_20.mat')
% plot(P_S_dB,EC_ERA_sim,'bo:'); hold on;
% plot(P_S_dB,EC_ERA_ana,'b-'); hold on;
% load('./data/data_Gamma_20_25_30.mat')
% plot(P_S_dB,EC_ERA_sim,'gs:'); hold on;
% plot(P_S_dB,EC_ERA_ana,'g-'); hold on;
% 
% % for OST
% 
% load('./data/data_LogNormal_25_25_25.mat')
% plot(P_S_dB,EC_ORA_LN_sim,'r^:'); hold on;
% plot(P_S_dB,EC_ORA_LN_ana,'r-'); hold on;
% load('./data/data_LogNormal_30_25_20.mat')
% plot(P_S_dB,EC_ORA_LN_sim,'bo:'); hold on;
% plot(P_S_dB,EC_ORA_LN_ana,'b-'); hold on;
% load('./data/data_LogNormal_20_25_30.mat')
% plot(P_S_dB,EC_ORA_LN_sim,'gs:'); hold on;
% plot(P_S_dB,EC_ORA_LN_ana,'g-'); hold on;
% 
% xlabel('Transmit SNR, $\bar{\rho}$ (dB)', 'Interpreter', 'Latex');
% ylabel('Ergodic capacity, $\bar{C}$ (b/s/Hz)', ...
%     'Interpreter', 'Latex');
% legend('$[25, 25, 25]$ (sim.)', '$[25, 25, 25]$ (ana.)',...
%     '$[30, 25. 20]$ (sim.)', '$[30, 25. 20]$ (ana.)',...
%     '$[20, 25. 30]$ (sim.)', '$[20. 25. 30]$ (ana.)',...
%     'Location','nw',...
%     'Interpreter', 'Latex');
% set(gca,'fontsize',16);

% figure;
% hold on; box on; grid on;
% plot(d,10*log10(beta_3GPP_LOS),'k-.','LineWidth',2);
% plot(d,10*log10(beta_3GPP_NLOS),'r--','LineWidth',2);
% xlabel('Distance $d$ [m]','Interpreter','Latex');
% ylabel('Channel gain $\beta(d)$ [dB]','Interpreter','Latex');
% legend({'UMi-LOS','UMi-NLOS'},'Interpreter','Latex');
% set(gca,'fontsize',18);
% ylim([-110 -50]);
