% Multi-RIS-aided Wireless Systems: Statistical Characterization and Performance Analysis
% Tri Nhu Do, Georges Kaddoum, Thanh Luan Nguyen, Daniel Benevides da Costa, and Zygmunt J. Haas
% https://arxiv.org/abs/2104.01912
% Version: 2021-04-05


%Script for Fig. 3

clear;
close all;

N_RIS = 5; %Number of RISs
L = 25*ones(1, N_RIS); %Element setting L1
P_S_dB = -5:25; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

%% Plot

load('data_Gamma_settingL1.mat')

subplot(1, 4, 1);

semilogy(P_S_dB, OP_ERA_Gamma_sim, 'r^'); hold on;
semilogy(P_S_dB, OP_ERA_Gamma_ana, 'r-'); hold on;

semilogy(P_S_dB, OP_ORA_Gamma_sim, 'bs'); hold on;
semilogy(P_S_dB, OP_ORA_Gamma_ana, 'b-'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Outage probability', 'Interpreter', 'Latex');
legend('ERA (sim.)', ...
    'ERA (ana. Eq.(xx))', ...
    'ORA (sim.)', ...
    'ORA (ana. Eq.(xx))', ...
    'Location', 'ne', ...
    'Interpreter', 'Latex');
axis([-5 25 10^(-5) 1e0]);
xticks([-5 5 15 25])
set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca, 'fontsize', 13);

%--------------------------------------------------------------------------
subplot(1, 4, 2);

plot(P_S_dB(1:3:30), EC_ERA_Gamma_sim(1:3:30), 'r^'); hold on;
plot(P_S_dB, EC_ERA_Gamma_ana, 'r-'); hold on;

plot(P_S_dB(1:3:30), EC_ORA_Gamma_sim(1:3:30), 'bs'); hold on;
plot(P_S_dB, EC_ORA_Gamma_ana, 'b-'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [b/s/Hz]', 'Interpreter', 'Latex');
legend('ERA (sim.)', ...
    'ERA (ana. Eq.(xx))', ...
    'ORA (sim.)', ...
    'ORA (ana. Eq.(xx))', ...
    'Location', 'nw', ...
    'Interpreter', 'Latex');
axis([-5 25 -Inf Inf]);
xticks([-5 5 15 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca, 'fontsize', 13);

%--------------------------------------------------------------------------
load('data_LogNormal_settingL1.mat')

subplot(1, 4, 3);

semilogy(P_S_dB, OP_ERA_LN_sim, 'r^'); hold on;
semilogy(P_S_dB, OP_ERA_LN_ana, 'r-'); hold on;

semilogy(P_S_dB, OP_ORA_LN_sim, 'bs'); hold on;
semilogy(P_S_dB, OP_ORA_LN_ana, 'b-'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Outage probability', 'Interpreter', 'Latex');
legend('ERA (sim.)', ...
    'ERA (ana. Eq.(xx))', ...
    'ORA (sim.)', ...
    'ORA (ana. Eq.(xx))', ...
    'Location', 'ne', ...
    'Interpreter', 'Latex');
axis([-5 25 10^(-5) 1e0]);
xticks([-5 5 15 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca, 'fontsize', 13);

%--------------------------------------------------------------------------

subplot(1, 4, 4);

plot(P_S_dB(1:3:30), EC_ERA_LN_sim(1:3:30), 'r^'); hold on;
plot(P_S_dB, EC_ERA_LN_ana, 'r-'); hold on;

plot(P_S_dB(1:3:30), EC_ORA_LN_sim(1:3:30), 'bs'); hold on;
plot(P_S_dB, EC_ORA_LN_ana, 'b-'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [b/s/Hz]', 'Interpreter', 'Latex');
legend('ERA (sim.)', ...
    'ERA (ana. Eq.(xx))', ...
    'ORA (sim.)', ...
    'ORA (ana. Eq.(xx))', ...
    'Location', 'nw', ...
    'Interpreter', 'Latex');
axis([-5 25 -Inf Inf]);
xticks([-5 5 15 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca, 'fontsize', 13);



