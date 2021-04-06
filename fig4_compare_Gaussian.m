% Multi-RIS-aided Wireless Systems: Statistical Characterization and Performance Analysis
% Tri Nhu Do, Georges Kaddoum, Thanh Luan Nguyen, Daniel Benevides da Costa, and Zygmunt J. Haas
% https://arxiv.org/abs/2104.01912
% Version: 2021-04-05

%% Simulation parameters

clear;
close all;

sim_trials = 1e6; %Number of simulation trails

R_th = 1; %Predefined target spectral efficiency [b/s/Hz]
SNR_th = 2^R_th-1; %Predefined SNR threshold

N_RIS = 5; %Number of RISs

% L = randi([20, 30], 1, N_RIS); %Number of elements at each RIS
L = 25*ones(1, N_RIS); %Element setting L1

kappa_nl = 1; %Amplitude reflection coefficient

%Nakagami m parameter
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

% %Random location setting
% x_RIS = x_area_min + (x_area_max-x_area_min)*rand(N_RIS, 1); % [num_RIS x 1] vector
% y_RIS = y_area_min + (y_area_max-y_area_min)*rand(N_RIS, 1);

%Loacation setting D1
x_RIS = [7; 13; 41; 75; 93];
y_RIS = [2; 6; 8; 4; 3];


% x_RIS = [10; 37; 57; 69; 87]; % good for fig 2; not good for fig 7
% y_RIS = [6; 3; 5; 7; 8];

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

P_S_dB = -5:25; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

SNRdB = P_S_dB - sigma2dBm; %Average transmit SNR, dB = dBm - dBm

%% SIMULATION

%Direct channel h_0
Omg_0 = pathloss_NLOS(d_sd)*antenna_gain_S; %Omega of S->D link

h_0 = random('Naka', m_0, Omg_0, [1, sim_trials]);

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

%% Fitting ERA

pd_Gamma = fitdist(Z_ERA', 'Gamma');
a_Gamma = pd_Gamma.a;
b_Gamma = pd_Gamma.b;
data_Gamma = gamrnd(a_Gamma, b_Gamma, [1, sim_trials]);

pd_LogNormal = fitdist(Z_ERA', 'Lognormal'); % Z^2 fit to LogNormal
mu_LogNormal = pd_LogNormal.mu;
sigma_LogNormal = pd_LogNormal.sigma;
data_LogNormal = lognrnd(mu_LogNormal, sigma_LogNormal, [1, sim_trials]);

pd_Normal = fitdist(Z_ERA', 'Normal');
mu_Normal = pd_Normal.mu;
sigma_Normal = pd_Normal.sigma;
data_Normal = normrnd(mu_Normal, sigma_Normal, [1, sim_trials]);

% pd_Burr = fitdist(Z_ERA', 'Burr');
% alpha_Burr = pd_Burr.alpha;
% c_Burr = pd_Burr.c;
% k_Burr = pd_Burr.k;
% data_Burr = random('burr', alpha_Burr, c_Burr, k_Burr, 1, sim_trials);

% pd_Weibull = fitdist(Z_ERA', 'Weibull');
% A_Weibull = pd_Weibull.A;
% B_Weibull = pd_Weibull.B;
% data_Weibull = wblrnd(A_Weibull, B_Weibull, [1, sim_trials]);

%% Fitting ORA

pd_Gamma = fitdist(R_ORA', 'Gamma');
a_Gamma = pd_Gamma.a;
b_Gamma = pd_Gamma.b;
data_Gamma_ORA = gamrnd(a_Gamma, b_Gamma, [1, sim_trials]);

pd_LogNormal = fitdist(R_ORA', 'Lognormal'); % Z^2 fit to LogNormal
mu_LogNormal = pd_LogNormal.mu;
sigma_LogNormal = pd_LogNormal.sigma;
data_LogNormal_ORA = lognrnd(mu_LogNormal, sigma_LogNormal, [1, sim_trials]);

pd_Normal = fitdist(R_ORA', 'Normal');
mu_Normal = pd_Normal.mu;
sigma_Normal = pd_Normal.sigma;
data_Normal_ORA = normrnd(mu_Normal, sigma_Normal, [1, sim_trials]);

%% CDF 

% figure;
% 
% [y, x] = ecdf(Z_ERA);
% plot(x, y, 'r'); hold on
% 
% [y, x] = ecdf(data_Gamma);
% plot(x, y); hold on
% 
% [y, x] = ecdf(data_LogNormal);
% plot(x, y); hold on
% 
% [y, x] = ecdf(data_Normal);
% plot(x, y); hold on
% 
% [y, x] = ecdf(data_Burr);
% plot(x, y); hold on
% 
% [y, x] = ecdf(data_Weibull);
% plot(x, y); hold on

%% PDF 

% %--------------------------------------------------------------------------
% %ORA
% 
% % figure;
% 
% [f_true, x_true] = ksdensity(R_ORA);
% plot(x_true, f_true, 'k--', 'linewidth', 2); 
% hold on
% 
% [f_Gamma, x_Gamma] = ksdensity(data_Gamma_ORA);
% plot(x_Gamma, f_Gamma, 'r--', 'linewidth', 2); 
% hold on
% 
% [f_LN, x_LN] = ksdensity(data_LogNormal_ORA);
% plot(x_LN, f_LN, 'g--', 'linewidth', 2);
% hold on
% 
% [f_Normal, x_Normal] = ksdensity(data_Normal_ORA);
% plot(x_Normal, f_Normal, 'b--', 'linewidth', 2);
% hold on

%--------------------------------------------------------------------------

%ERA

subplot(1, 3, 1);

[f_true, x_true] = ksdensity(Z_ERA);
plot(x_true, f_true, 'k-', 'linewidth', 2); 
hold on

[f_Gamma, x_Gamma] = ksdensity(data_Gamma);
plot(x_Gamma, f_Gamma, 'r-', 'linewidth', 2); 
hold on

[f_LN, x_LN] = ksdensity(data_LogNormal);
plot(x_LN, f_LN, 'g-', 'linewidth', 2);
hold on

[f_Normal, x_Normal] = ksdensity(data_Normal);
plot(x_Normal, f_Normal, 'b-', 'linewidth', 2);
hold on


xlabel('$x$', 'Interpreter', 'Latex');
ylabel('PDF', 'Interpreter', 'Latex');
legend('True (sim.)', ...
    'Gamma (ana.)', ...
    'LogNormal (ana.)', ...
    'Gaussian (sim.)', ...
    'Location', 'se',...
    'Interpreter', 'Latex');
% axis([0.5*1e-5, 2*1e-5, -Inf, Inf])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

zoomPlot(x_true,f_true,[2.0e-5, 2.2e-5],[1e-5 .1 0.3 0.2],[0 0]); hold on
plot(x_true, f_true, 'k-', 'linewidth', 2); hold on
plot(x_Gamma, f_Gamma, 'r-', 'linewidth', 2); hold on
plot(x_LN, f_LN, 'g-', 'linewidth', 2);hold on
plot(x_Normal, f_Normal, 'b-', 'linewidth', 2);hold on

%% OP

for isnr = 1:length(SNRdB)
    snr = 10^(SNRdB(isnr)/10);
    %
    OP_ERA_sim(isnr) = mean(snr*Z2_ERA < SNR_th);
    OP_ERA_Gamma(isnr) = mean(snr*(data_Gamma.^2) < SNR_th);
    OP_ERA_LN(isnr) = mean(snr*(data_LogNormal.^2) < SNR_th);
    OP_ERA_Normal(isnr) = mean(snr*(data_Normal.^2) < SNR_th);
    
    OP_ORA_sim(isnr) = mean(snr*R2_ORA < SNR_th);
    OP_ORA_Gamma(isnr) = mean(snr*(data_Gamma_ORA.^2) < SNR_th);
    OP_ORA_LN(isnr) = mean(snr*(data_LogNormal_ORA.^2) < SNR_th);
    OP_ORA_Normal(isnr) = mean(snr*(data_Normal_ORA.^2) < SNR_th);

    fprintf('Outage Probability, SNR = %d \n', SNRdB(isnr));
end

load('./data/data_Gamma_settingL1.mat')
load('./data/data_LogNormal_settingL1.mat')

subplot(1, 3, 2);

semilogy(P_S_dB, OP_ERA_sim, 'k-', 'linewidth', 2); hold on;
semilogy(P_S_dB, OP_ERA_Gamma_ana, 'r-', 'linewidth', 2); hold on;
semilogy(P_S_dB, OP_ERA_LN_ana, 'g-', 'linewidth', 2); hold on;
semilogy(P_S_dB, OP_ERA_Normal, 'b-', 'linewidth', 2); hold on;

semilogy(P_S_dB, OP_ORA_sim, 'k-', 'linewidth', 2); hold on;
semilogy(P_S_dB, OP_ORA_Gamma_ana, 'r-', 'linewidth', 2); hold on;
semilogy(P_S_dB, OP_ORA_LN_ana, 'g-', 'linewidth', 2); hold on;
semilogy(P_S_dB, OP_ORA_Normal, 'b-', 'linewidth', 2); hold on;


% zoomPlot_semilogy(P_S_dB,OP_ERA_sim,[1, 2],[10 .1 0.4 0.3],[0 0]); hold on
% semilogy(P_S_dB, OP_ERA_sim, 'k-', 'linewidth', 2); hold on;
% semilogy(P_S_dB, OP_ERA_Gamma_ana, 'r-', 'linewidth', 2); hold on;
% semilogy(P_S_dB, OP_ERA_LN_ana, 'g-', 'linewidth', 2); hold on;
% semilogy(P_S_dB, OP_ERA_Normal, 'b-', 'linewidth', 2); hold on;


xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Outage probability', 'Interpreter', 'Latex');
legend('True (sim.)', ...
    'Gamma (ana.)', ...
    'LogNormal (ana.))', ...
    'Gaussian (sim.)', ...
    'Location', 'ne',...
    'Interpreter', 'Latex');
axis([-Inf Inf 10^(-5) 10^(0)]);
xticks([-5 5 15 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);



%% EC

for isnr = 1:length(SNRdB)
    snr = 10^(SNRdB(isnr)/10);

    EC_ERA_sim(isnr) = mean(log2(1+snr*Z2_ERA));
    EC_ERA_Gamma(isnr) = mean(log2(1+snr*(data_Gamma.^2)));
    EC_ERA_LN(isnr) = mean(log2(1+snr*(data_LogNormal.^2)));
    EC_ERA_Normal(isnr) = mean(log2(1+snr*(data_Normal.^2)));
    
    EC_ORA_sim(isnr) = mean(log2(1+snr*R2_ORA));
    EC_ORA_Gamma(isnr) = mean(log2(1+snr*(data_Gamma_ORA.^2)));
    EC_ORA_LN(isnr) = mean(log2(1+snr*(data_LogNormal_ORA.^2)));
    EC_ORA_Normal(isnr) = mean(log2(1+snr*(data_Normal_ORA.^2)));

    fprintf('EC, SNR = %d \n', SNRdB(isnr));
end

% figure;
subplot(1, 3, 3);

plot(P_S_dB(1:3:31), EC_ERA_sim(1:3:31), 'k-', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ERA_Gamma_ana(1:3:31), 'ro:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ERA_LN_ana(1:3:31), 'g*:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ERA_Normal(1:3:31), 'b>:', 'linewidth', 1); hold on;

plot(P_S_dB(1:3:31), EC_ORA_sim(1:3:31), 'k-', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ORA_Gamma_ana(1:3:31), 'ro:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ORA_LN_ana(1:3:31), 'g*:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ORA_Normal(1:3:31), 'b>:', 'linewidth', 1); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [b/s/Hz]', 'Interpreter', 'Latex');
legend('True (sim.)', ...
    'Gamma (ana.)', ...
    'LogNormal (ana.)', ...
    'Gaussian (sim.)', ...
    'Location', 'ne',...
    'Interpreter', 'Latex');
axis([-5 25 -Inf Inf]);
xticks([-5 5 15 25])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

zoomPlot(P_S_dB,EC_ERA_sim,[13, 15],[.5 .1 0.4 0.3],[0 0]); hold on
plot(P_S_dB(1:3:31), EC_ERA_sim(1:3:31), 'k-', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ERA_Gamma_ana(1:3:31), 'ro:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ERA_LN_ana(1:3:31), 'g*:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ERA_Normal(1:3:31), 'b>:', 'linewidth', 1); hold on;

plot(P_S_dB(1:3:31), EC_ORA_sim(1:3:31), 'k-', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ORA_Gamma_ana(1:3:31), 'ro:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ORA_LN_ana(1:3:31), 'g*:', 'linewidth', 1); hold on;
plot(P_S_dB(1:3:31), EC_ORA_Normal(1:3:31), 'b>:', 'linewidth', 1); hold on;
