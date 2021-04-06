% Multi-RIS-aided Wireless Systems: Statistical Characterization and Performance Analysis
% Tri Nhu Do, Georges Kaddoum, Thanh Luan Nguyen, Daniel Benevides da Costa, and Zygmunt J. Haas
% https://arxiv.org/abs/2104.01912
% Version: 2021-04-05


clear;
close all;

sim_trials = 1e6; %Number of simulation trails

R_th = 3; %Predefined target SE [bit/s/Hz]
SNR_th = 2^R_th-1; %Predefined SNR threshold

N_RIS = 5; %Number of RISs

% L = randi([20, 30], 1, N_RIS); %Number of elements at each RIS
% L = 25*ones(1, N_RIS); %Element setting L1
% L = 40*ones(1, N_RIS); %Element setting L2
% L = [20 30 40 50 60]; %Element setting L3
L = [60 50 40 30 20]; %Element setting L4

kappa_nl = 1; %Amplitude reflection coefficient

%Nakagami m parameter (random)
% m_0 = 2.5 + rand; %Scale parameter, Heuristic setting
m_0 = 3.551542983398091;

% m_h = 2.5 + rand(N_RIS, 1);
m_h = [2.72692063869142;3.18539808499718;2.97177594571026;2.97218021688479;2.98236942594706];

% m_g = 2.5 + rand(N_RIS, 1);
m_g = [3.24995839689103;2.63233191045873;3.16576836183012;3.40081750859119;3.45261429002219];

%--------------------------------------------------------------------------

%Network area
x_area_min = 0;
x_area_max = 100;
y_area_min = 0;
y_area_max = 10;

%Source location
x_source = x_area_min;
y_source = y_area_min;

%Destination location
x_des = x_area_max;
y_des = y_area_min;

%Random topology
x_RIS = x_area_min + (x_area_max-x_area_min)*rand(N_RIS, 1); % [num_RIS x 1] vector
y_RIS = y_area_min + (y_area_max-y_area_min)*rand(N_RIS, 1);

%Location setting D1
x_RIS = [7; 13; 41; 75; 93];
y_RIS = [2; 6; 8; 4; 3];

%Compute location of nodes
pos_source = [x_source, y_source];
pos_des = [x_des, y_des];
pos_RIS = [x_RIS, y_RIS]; %[num_RIS x 2] matrix

%Compute distances
d_sr = sqrt(sum((pos_source - pos_RIS).^2 , 2)); %[num_RIS x 1] vector
d_rd = sqrt(sum((pos_RIS - pos_des).^2 , 2));
d_sd = sqrt(sum((pos_source - pos_des).^2 , 2));

%--------------------------------------------------------------------------
%Path-loss model
%Carrier frequency (in GHz)
fc = 3; % GHz

%3GPP Urban Micro in 3GPP TS 36.814, Mar. 2010. 
%Note that x is measured in meter

%NLoS path-loss component based on distance
pathloss_NLOS = @(x) db2pow(-22.7 - 26*log10(fc) - 36.7*log10(x));

antenna_gain_S = db2pow(5); %Source antenna gain, dBi
antenna_gain_RIS = db2pow(5); %Gain of each element of a RIS, dBi
antenna_gain_D = db2pow(0); %Destination antenna gain, dBi

%--------------------------------------------------------------------------
%Noise power and Transmit power P_S
%Bandwidth
BW = 10e6; %10 MHz

%Noise figure (in dB)
noiseFiguredB = 10;

%Compute the noise power in dBm
sigma2dBm = -174 + 10*log10(BW) + noiseFiguredB; %-94 dBm
sigma2 = db2pow(sigma2dBm);

P_S_dB = -5:.3:25; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

SNRdB = P_S_dB - sigma2dBm; %Average transmit SNR, dB = dBm - dBm

%% SIMULATION

%Direct channel h_0
Omg_0 = pathloss_NLOS(d_sd)*antenna_gain_S; %Omega of S->D link

h_0 = random('Naka', m_0, Omg_0, [1, sim_trials]);
SNR_h0 = abs(h_0).^2;

V_n = zeros(N_RIS, sim_trials);
Omg_h = zeros(N_RIS, 1); 
Omg_g = zeros(N_RIS, 1); 

for nn = 1:N_RIS
    for kk = 1:L(nn)
        Omg_h(nn) = pathloss_NLOS(d_sr(nn))*antenna_gain_S*antenna_gain_RIS*L(nn); %Omega S->R
        Omg_g(nn) = pathloss_NLOS(d_rd(nn))*antenna_gain_RIS*L(nn)*antenna_gain_D; %Omega R->D
        
        h_nl = random('Naka', m_h(nn), Omg_h(nn), [1, sim_trials]);
        g_nl = random('Naka', m_g(nn), Omg_g(nn), [1, sim_trials]);
        
        U_nl = kappa_nl * h_nl .* g_nl;
        
        V_n(nn, :) = V_n(nn, :) + U_nl;
    end
end

%ERA scheme
T_ERA   = sum(V_n, 1);
Z_ERA   = h_0 + T_ERA; %Magnitude of the e2e channel
Z2_ERA  = Z_ERA.^2; %Squared magnitude of the e2e channel

%ORA scheme
V_M_ORA = max(V_n, [], 1); %V_M for the best RIS
R_ORA   = h_0 + V_M_ORA; %Magnitude of the e2e channel
R2_ORA  = R_ORA.^2; %Squared magnitude of the e2e channel

%% OUTAGE PROBABILITY

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); %i.e., 10^(SNRdB/10)
    
    OP_SISO(idx) = mean(avgSNR*SNR_h0 < SNR_th);
    OP_ERA_Gamma_sim(idx) = mean(avgSNR*Z2_ERA < SNR_th);
    OP_ORA_Gamma_sim(idx) = mean(avgSNR*R2_ORA < SNR_th);
    
    fprintf('Gamma, Outage probability, SNR = %d \n', round(SNRdB(idx)));
end

figure;
semilogy(P_S_dB, OP_SISO, 'gd:'); hold on;
semilogy(P_S_dB, OP_ERA_Gamma_sim, 'r*:'); hold on;
semilogy(P_S_dB, OP_ORA_Gamma_sim, 'bo:'); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Outage Probability, $P_{\rm out}$', 'Interpreter', 'Latex');
legend('SISO (sim.)',...
    'ERA (sim.)', ...
    'ORA (sim.)', ...
    'Location','SW',...
    'Interpreter', 'Latex');
axis([-Inf Inf 10^(-5) 10^(0)]);

%% ERGODIC CAPACITY

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); %10^(SNRdB(idx)/10)
    
    EC_SISO(idx) = mean(log2(1 + avgSNR*SNR_h0));
    EC_ERA_Gamma_sim(idx) = mean(log2(1+avgSNR*Z2_ERA));
    EC_ORA_Gamma_sim(idx) = mean(log2(1+avgSNR*R2_ORA));
    
    fprintf('Gamma, Ergodic capacity, SNR = %d \n', round(SNRdB(idx)));
end

figure;
plot(P_S_dB, EC_SISO, 'go:'); hold on;
plot(P_S_dB, EC_ERA_Gamma_sim, 'ro:'); hold on;
plot(P_S_dB, EC_ORA_Gamma_sim, 'bs:'); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Ergodic Capacity (bit/s/Hz)', 'Interpreter', 'Latex');
legend('SISO (sim.)',...
    'ERA (sim.)',...
    'ORA (sim.)',...
    'Interpreter', 'Latex',...
    'Location','NW');

%% SAVE DATA

save('data_Gamma_setL4.mat',...
    'OP_SISO',...
    'OP_ERA_Gamma_sim',...
    'OP_ORA_Gamma_sim',...
    'EC_SISO',...
    'EC_ERA_Gamma_sim',...
    'EC_ORA_Gamma_sim')
