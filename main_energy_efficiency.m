clear all
close all
%
sim_times = 1e5; % number of simulation trails
%
snrdB = -20:.1:40; % average transmit SNR
r_th = 1.0; % predefined data rate bits/s/Hz
snr_th = 2^(r_th)-1; % predefined SNR threshold
%
N_RIS = 3; % number of RISs
L = [30 25 20]; % number of elements at each RIS
%
kappa_nl = 1;
% m_h = linspace(2, 3, N_RIS)';
% m_g = linspace(3, 2, N_RIS)';
m0 = 1.5;
m_h = [3.3, 2.9, 2.3]';
m_g = [2.1, 2.7, 3.5]';
%
loss = 1e3; % path-loss at reference distance d0
epsilon = 4; % path-loss exponent
d0 = 1; % reference distance in meter
%
% network area
x_area_min = 0;
x_area_max = 500;
y_area_min = 0;
y_area_max = 20;
%
% source location
x_source = x_area_min;
y_source = y_area_min;
%
% destination location
x_des = x_area_max;
y_des = y_area_min;
%
% random locations of relays
% x_ris = x_area_min + (x_area_max-x_area_min)*rand(N_RIS,1); % [num_relay x 1] vector
% y_ris = y_area_min + (y_area_max-y_area_min)*rand(N_RIS,1);

% x_ris = x_area_min + (x_area_max-x_area_min)*linspace(0.25, 0.75, N_RIS)'; % [num_relay x 1] vector
% y_ris = .7* ones(N_RIS,1)* y_area_max;

x_ris = ones(N_RIS,1); x_ris(1) = 70; x_ris(2) = 230; x_ris(3) = 450;
y_ris = ones(N_RIS,1); y_ris(1) = 17; y_ris(2) = 5; y_ris(3) = 10;

% location
pos_source = [x_source, y_source];
pos_des = [x_des, y_des];
pos_relay = [x_ris, y_ris]; % [num_relay x 2] matrix

% distances
d_sr = sqrt(sum((pos_source - pos_relay).^2 , 2)); % [num_relay x 1] vector
d_rd = sqrt(sum((pos_relay - pos_des).^2 , 2));
d_sd = sqrt(sum((pos_source - pos_des).^2 , 2));

% non-identical Omega
Omgh = loss./((d_sr./d0).^(epsilon/2)); % [num_relay x 1] vector
Omgg = loss./((d_rd./d0).^(epsilon/2));
Omg0 = loss./((d_sd./d0).^(epsilon/2));

%% Simulation
h_0 = random('Naka',m0,Omg0,[1,sim_times]);

V_n = zeros(N_RIS, sim_times);
for n = 1:N_RIS
    for l = 1:L(n)
        h_nl = random('Naka', m_h(n), Omgh(n), [1,sim_times]);
        g_nl = random('Naka', m_g(n), Omgg(n), [1,sim_times]);
        U_nl = kappa_nl * h_nl .* g_nl;
        %
        V_n(n,:) = V_n(n,:) + U_nl;
    end
end
%
T_cct   = sum(V_n,1);
Z_cct   = h_0 + T_cct;
Z2_cct  = Z_cct.^2;
%
V_M_ost = max(V_n,[],1); % V_M : best RIS
R_ost   = h_0 + V_M_ost;
R2_ost  = R_ost.^2;


%% Performance metrics

for isnr = 1:length(snrdB)
    fprintf('OP, SNR = %d \n', snrdB(isnr));
    snr = 10^(snrdB(isnr)/10);
    
    EC_CCT_sim(isnr) = mean(log2(1+snr*Z2_cct));
    
    EC_OST_sim(isnr) = mean(log2(1+snr*R2_ost));
end

%% Numerical Performance Analysis

% Ergodic Capacity

% figure;
% plot(snrdB, EC_CCT_sim, 'ro-'); hold on;
% plot(snrdB, EC_OST_sim, 'b+-'); hold on;
% xlabel('Transmit SNR (dB))');
% ylabel('EC (b/s//Hz');
% legend('CCT (sim.)',...
%     'OST (sim.)',...
%     'Location','Best');

% Transmit Power
EC = 1:15;
P_tx_CCT = zeros(1, length(EC));

P_tx_OST = zeros(1, length(EC));
for ii = 1:length(EC)
    idx = find(EC_CCT_sim <= EC(ii));
    P_tx_CCT(ii) = 10^(snrdB(idx(end))/10);
    idx = find(EC_OST_sim <= EC(ii));
    P_tx_OST(ii) = 10^(snrdB(idx(end))/10);
end

figure;
plot(EC, P_tx_CCT, 'ro-'); hold on
plot(EC, P_tx_OST, 'bs-'); hold on

xlabel('Achievable rate (b/s/Hz)')
ylabel('Transmit power (mW)')
legend('CCT', 'OST', 'location', 'best')

% Engery Efficiency

P_c_S = 10^(10/10); % Circuit dissipated power at S = 10 dBm
P_c_D = 10^(10/10); %% Circuit dissipated power at S = 10 dBm
P_c_element = 7.8; % mW
BW = 180*1e3; % 180kHz

for kk=1:length(EC)
    P_total_CCT(kk) = P_tx_CCT(kk) + 75*P_c_element + P_c_S + P_c_D;
    P_total_OST(kk) = P_tx_OST(kk) + 25*P_c_element + P_c_S + P_c_D;
end

figure
plot(EC, BW*EC./P_total_CCT/1e3, 'ro-'); hold on
plot(EC, BW*EC./P_total_OST/1e3, 'b*-'); hold on
xlabel('Achievable rate (bit/s/Hz)')
ylabel('Energy efficiency (Mbit/Joule)')
legend('CCT', 'OST', 'location', 'nw')

% x0 = 100; y0 = 100; width = 400; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====

% %% Network Topology
% 
% %-----Topo
% fig_topo = figure;
% scatter(x_source, y_source, 'bo', 'filled')
% hold on
% scatter(x_des, y_des, 'go', 'filled')
% hold on
% scatter(x_ris, y_ris, 'rs', 'filled')
% hold on
% 
% for kk = 1:N_RIS
%     text(x_ris(kk), y_ris(kk), num2str(kk));
%     hold on
% end
% 
% axis([x_area_min x_area_max y_area_min y_area_max])
% legend('Source', 'Destination', 'RISs')
% 
% x0 = 100; y0 = 100; width = 400; height = 200;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
% %=====DONE
% %