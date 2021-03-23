clear;
close all;

sim_trials = 1e6; %Number of simulation trails
R_th = 1; %Predefined data rate bit/s/Hz
SNR_th = 2^R_th-1; %Predefined SNR threshold

N_RIS = 5; %Number of RISs
L = 25*ones(1, N_RIS); %Number of elements at each RIS
% L = randi([20, 30], 1, N_RIS);

kappa_nl = 1; %Amplitude reflection coefficient

%Nakagami m parameter (random)
m_0 = 2.5 + rand; %Scale parameter, Heuristic setting
% m_0 = 3.551542983398091;

m_h = 2.5 + rand(N_RIS, 1);
% m_h = [3.566804625889170;3.588550582888764;3.364073718400015];

m_g = 2.5 + rand(N_RIS, 1);
% m_g = [3.526226429649985; 3.564715591043445; 3.502113633564524];
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

% x_RIS = [13; 52; 79];
x_RIS = [10; 37; 57; 69; 87];
% y_RIS = [7; 5; 8];
y_RIS = [6; 3; 5; 7; 8];

% %Line topology
% x_RIS = linspace(x_area_min+x_area_max/N_RIS, x_area_max-x_area_max/N_RIS, N_RIS)'; % [num_RIS x 1] vector
% y_RIS = y_area_max*ones(N_RIS, 1);

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

P_S_dB = 0:30; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

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

%ERA

% figure;
subplot(1, 3, 1);


[f_true, x_true] = ksdensity(Z_ERA);
plot(x_true, f_true, 'k-', 'linewidth', 2); 
hold on

[f_Gamma, x_Gamma] = ksdensity(data_Gamma);
plot(x_Gamma, f_Gamma, 'r-', 'linewidth', 2); 
hold on

[f_Lognormal, x_Lognormal] = ksdensity(data_LogNormal);
plot(x_Lognormal, f_Lognormal, 'g-', 'linewidth', 2);
hold on

[f_Normal, x_Normal] = ksdensity(data_Normal);
plot(x_Normal, f_Normal, 'b-', 'linewidth', 2);
hold on

%--------------------------------------------------------------------------
%ORA

% figure;

[f_true, x_true] = ksdensity(R_ORA);
plot(x_true, f_true, 'k--', 'linewidth', 2); 
hold on

[f_Gamma, x_Gamma] = ksdensity(data_Gamma_ORA);
plot(x_Gamma, f_Gamma, 'r--', 'linewidth', 2); 
hold on

[f_Lognormal, x_Lognormal] = ksdensity(data_LogNormal_ORA);
plot(x_Lognormal, f_Lognormal, 'g--', 'linewidth', 2);
hold on

[f_Normal, x_Normal] = ksdensity(data_Normal_ORA);
plot(x_Normal, f_Normal, 'b--', 'linewidth', 2);
hold on

xlabel('$x$', 'Interpreter', 'Latex');
ylabel('PDF', 'Interpreter', 'Latex');
legend('True (sim.)', ...
    'Gamma (ana. Thm. xx)', ...
    'LogNormal (ana. Thm. xx)', ...
    'Gaussian (sim.)', ...
    'Location', 'se',...
    'Interpreter', 'Latex');
axis([0.5*1e-5, 2*1e-5, -Inf, Inf])

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);


%% OP

for isnr = 1:length(SNRdB)
    snr = 10^(SNRdB(isnr)/10);
    %
    Pout_ERA_sim(isnr) = mean(snr*Z2_ERA < SNR_th);
    Pout_ERA_Gamma(isnr) = mean(snr*(data_Gamma.^2) < SNR_th);
    Pout_ERA_Lognormal(isnr) = mean(snr*(data_LogNormal.^2) < SNR_th);
    Pout_ERA_Normal(isnr) = mean(snr*(data_Normal.^2) < SNR_th);
    
    Pout_ORA_sim(isnr) = mean(snr*R2_ORA < SNR_th);
    Pout_ORA_Gamma(isnr) = mean(snr*(data_Gamma_ORA.^2) < SNR_th);
    Pout_ORA_Lognormal(isnr) = mean(snr*(data_LogNormal_ORA.^2) < SNR_th);
    Pout_ORA_Normal(isnr) = mean(snr*(data_Normal_ORA.^2) < SNR_th);

    fprintf('Outage Probability, SNR = %d \n', SNRdB(isnr));
end

% figure;
subplot(1, 3, 2);

semilogy(P_S_dB, Pout_ERA_sim, 'k-', 'linewidth', 2); hold on;
semilogy(P_S_dB, Pout_ERA_Gamma, 'r-', 'linewidth', 2); hold on;
semilogy(P_S_dB, Pout_ERA_Lognormal, 'g-', 'linewidth', 2); hold on;
semilogy(P_S_dB, Pout_ERA_Normal, 'b-', 'linewidth', 2); hold on;


semilogy(P_S_dB, Pout_ORA_sim, 'k--', 'linewidth', 2); hold on;
semilogy(P_S_dB, Pout_ORA_Gamma, 'r--', 'linewidth', 2); hold on;
semilogy(P_S_dB, Pout_ORA_Lognormal, 'g--', 'linewidth', 2); hold on;
semilogy(P_S_dB, Pout_ORA_Normal, 'b--', 'linewidth', 2); hold on;

xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Outage Probability', 'Interpreter', 'Latex');
legend('True (sim.)', ...
    'Gamma (ana. Eq. (xx))', ...
    'LogNormal (ana. Eq. (xx))', ...
    'Gaussian (sim.)', ...
    'Location', 'se',...
    'Interpreter', 'Latex');
axis([-Inf Inf 10^(-5) 10^(0)]);

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%% EC

% EC_ERA_Gamma_sim(idx) = mean(log2(1+avgSNR*Z2_ERA));


for isnr = 1:length(SNRdB)
    snr = 10^(SNRdB(isnr)/10);

    EC_ERA_sim(isnr) = mean(log2(1+snr*Z2_ERA));
    EC_ERA_Gamma(isnr) = mean(log2(1+snr*(data_Gamma.^2)));
    EC_ERA_Lognormal(isnr) = mean(log2(1+snr*(data_LogNormal.^2)));
    EC_ERA_Normal(isnr) = mean(log2(1+snr*(data_Normal.^2)));
    
    EC_ORA_sim(isnr) = mean(log2(1+snr*R2_ORA));
    EC_ORA_Gamma(isnr) = mean(log2(1+snr*(data_Gamma_ORA.^2)));
    EC_ORA_Lognormal(isnr) = mean(log2(1+snr*(data_LogNormal_ORA.^2)));
    EC_ORA_Normal(isnr) = mean(log2(1+snr*(data_Normal_ORA.^2)));

    fprintf('EC, SNR = %d \n', SNRdB(isnr));
end

% figure;
subplot(1, 3, 3);

plot(P_S_dB, EC_ERA_sim, 'k-', 'linewidth', 2); hold on;
plot(P_S_dB, EC_ERA_Gamma, 'ro:', 'linewidth', 2); hold on;
plot(P_S_dB, EC_ERA_Lognormal, 'g*:', 'linewidth', 2); hold on;
plot(P_S_dB, EC_ERA_Normal, 'b>:', 'linewidth', 2); hold on;

plot(P_S_dB, EC_ORA_sim, 'k-', 'linewidth', 2); hold on;
plot(P_S_dB, EC_ORA_Gamma, 'ro:', 'linewidth', 2); hold on;
plot(P_S_dB, EC_ORA_Lognormal, 'g*:', 'linewidth', 2); hold on;
plot(P_S_dB, EC_ORA_Normal, 'b>:', 'linewidth', 2); hold on;



xlabel('$P_{\rm S}$ (dBm)', 'Interpreter', 'Latex');
ylabel('Ergodic Capacity (bit/sec/Hz)', 'Interpreter', 'Latex');
legend('True (sim.)', ...
    'Gamma (ana. Eq. (xx))', ...
    'LogNormal (ana. Eq. (xx))', ...
    'Gaussian (sim.)', ...
    'Location', 'se',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);


set(gca,'LooseInset',get(gca,'TightInset'))