clear;
close all;

sim_trials = 1e5; %Number of simulation trails


N_RIS = 5; %Number of RISs

% L = randi([20, 30], 1, N_RIS); %Number of elements at each RIS
% L = 25*ones(1, N_RIS); %Element setting L1
% L = 40*ones(1, N_RIS); %Element setting L2
% L = [20 30 40 50 60]; %Element setting L3
L = [60 50 40 30 20]; %Element setting L4

kappa_nl = 1; %Amplitude reflection coefficient


%Nakagami m parameter (random)
% m_0 = 2.5 + rand; %Scale parameter, Heuristic setting
m_0 = 3.236971152246479;

% m_h = 2.5 + rand(N_RIS, 1);
m_h = [3.40155918408956;3.15391856405791;3.29185084374551;2.93320811970437;3.21687693990103;2.56560300723790;3.07496768345179;3.33695541612840;3.00097856649605;3.11363316470094;2.89279112416316;3.07768794734794;3.02478972515711;2.89111670381469;2.67843264888998;3.17110899438754;3.04047559088993;3.28466155980647;3.48297567696655;3.03207250187595];

% m_g = 2.5 + rand(N_RIS, 1);
m_g = [2.58039894187387;2.64224899594153;3.11510661856003;3.38410720098219;3.39574370452735;2.65671648677722;3.40919813750129;3.42219275398468;2.88954481804285;3.45714034480911;3.40478234028093;3.43050220356429;2.75184550287702;3.06344191805314;2.69755997249909;2.63740256620508;3.04720955794659;2.67421093631393;2.89889509144691;2.87846927541581];

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

% %Location setting D1
% x_RIS = [7; 13; 41; 75; 93];
% y_RIS = [2; 6; 8; 4; 3];

%Location setting D2
x_RIS = [5; 13; 37; 69; 91];
y_RIS = [2; 7; 6; 1; 3];

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

P_S_dB = 23; %Transmit power of the source, dBm, e.g., 200mW = 23dBm

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

%% ERGODIC CAPACITY

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); %10^(SNRdB(idx)/10)
    
    EC_ERA_Gamma_sim(idx) = mean(log2(1+avgSNR*Z2_ERA));
    
end

% figure;
subplot(2,1,2);

plot(P_S_dB, EC_ERA_Gamma_sim, 'ro:'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [bit/s/Hz]', 'Interpreter', 'Latex');
legend('ERA (sim.)',...
    'Interpreter', 'Latex',...
    'Location','NW');

%% Network topology

% figure;
subplot(2,1,1);

scatter(x_source, y_source, 100, 'bo', 'filled')
hold on
scatter(x_des, y_des, 100, 'gd', 'filled')
hold on
scatter(x_RIS, y_RIS, 100, 'rs', 'filled')
hold on

for kk = 1:N_RIS
    text(x_RIS(kk)+3, y_RIS(kk), num2str(kk));
    hold on
end

xlabel('$d_{\rm SD}$ [m]', 'Interpreter', 'Latex')
ylabel('$H$ [m]', 'Interpreter', 'Latex')
axis([x_area_min x_area_max y_area_min y_area_max])
legend('$\rm S$', '$\rm D$', '$\mathrm{R}_i$',...
    'Interpreter', 'Latex',...
    'Location', 'ne')

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

disp(EC_ERA_Gamma_sim)