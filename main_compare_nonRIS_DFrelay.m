clear all
close all
%
sim_times = 1e6; % number of simulation trails
%
snrdB = 0:2:40; % average transmit SNR
R_th = 1.0; % predefined data rate bits/s/Hz
snr_th = 2^(R_th)-1; % predefined SNR threshold

%
N_RIS = 3; % number of RISs
L = [10 10 10]; % number of elements at each RIS
%

snr_th_HD = 2^(2*R_th) - 1;
snr_th_FD = 2^R_th - 1;

kappa_nl = 1;
% m_h = linspace(1, 3, N_RIS)';
% m_g = linspace(1, 3, N_RIS)';
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

% x_ris = x_area_min + (x_area_max-x_area_min)*linspace(.1, .9, N_RIS)'; % [num_relay x 1] vector
% y_ris = .7* ones(N_RIS,1)* y_area_max;

x_ris = ones(N_RIS,1); x_ris(1) = 70; x_ris(2) = 230; x_ris(3) = 450;
y_ris = ones(N_RIS,1); y_ris(1) = 17; y_ris(2) = 5; y_ris(3) = 10;

% location
pos_source = [x_source, y_source];
pos_des = [x_des, y_des];
pos_ris = [x_ris, y_ris]; % [num_relay x 2] matrix

% distances
d_sr = sqrt(sum((pos_source - pos_ris).^2 , 2)); % [num_relay x 1] vector
d_rd = sqrt(sum((pos_ris - pos_des).^2 , 2));
d_sd = sqrt(sum((pos_source - pos_des).^2 , 2));

% non-identical Omega
Omgh = loss./((d_sr./d0).^(epsilon/2)); % [num_relay x 1] vector
Omgg = loss./((d_rd./d0).^(epsilon/2));
Omg0 = loss./((d_sd./d0).^(epsilon/2));

%% Simulation
h_0 = random('Naka',m0,Omg0,[1,sim_times]);

V_n = zeros(N_RIS, sim_times);
h_SR = [];
h_RD = [];
for n = 1:N_RIS
    for l = 1:L(n)
        h_nl = random('Naka', m_h(n), Omgh(n), [1,sim_times]);
        g_nl = random('Naka', m_g(n), Omgg(n), [1,sim_times]);
        U_nl = kappa_nl * h_nl .* g_nl;
        %
        V_n(n,:) = V_n(n,:) + U_nl;
        h_SR = [h_SR; h_nl];
        h_RD = [h_RD; g_nl];
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

%
% For conventional dual-hop relaying
snr_h0 = (abs(h_0).^2); % for direction S-D link

h_SR = random('Naka', m_h(2), Omgh(2), [1,sim_times]);
h_RD = random('Naka', m_g(2), Omgg(2), [1,sim_times]);

% HD-DF relaying
snr_SR = abs(h_SR).^2;
snr_RD = abs(h_RD).^2;

% Full-Duplex (FD) relaying

% Rician parameters
K_dB        = 5; % Rician K-factor in dB
K           = 10.^(K_dB./10);
lrsi        = K/(K+1);                            % lambda_rsi
% Rician Distribution parameters
chi          = sqrt(K/(K+1)*lrsi); % Noncentrality parameter
sigma       = sqrt(lrsi/2/(K+1)); % Scale parameter

h_LI = random('rician',chi,sigma,[1,sim_times]);
snr_LI = abs(h_LI).^2;

%% Outage Probability

%-----OP
for isnr = 1:length(snrdB)
    snr = 10^(snrdB(isnr)/10);
    fprintf('SNR = %d dB\n', snrdB(isnr));
    
    % OP
    % CCT scheme
    P_out_CCT_Gamma_sim(isnr) = mean(snr*Z2_cct < snr_th);
    
    % OST scheme
    P_out_OST_Gamma_sim(isnr) = mean(snr*R2_ost < snr_th);
    
    % non_RIS
    Pout_nonRIS(isnr) = mean(snr*(abs(h_0).^2) < snr_th);
    
    
    % EC
    EC_CCT_sim(isnr) = mean(log2(1 + snr*Z2_cct));
    
    EC_OST_sim(isnr) = mean(log2(1 + snr*R2_ost));
    
    EC_nonRIS(isnr) = mean(log2(1 + snr*(abs(h_0).^2)));
    
    % HD-AF
    snr_e2e_HD_AF = snr*snr_h0 + (snr.*snr_SR.*snr.*snr_RD) ...
        ./(snr.*snr_SR + snr.*snr_RD + 1);
    
    Pout_HD_AF(isnr) = mean(snr_e2e_HD_AF < snr_th_HD);
    
    EC_HD_AF(isnr) = mean(log2(1 + snr_e2e_HD_AF))/2;
    
    % HD-DF
    snr_e2e_HD_DF =  snr*snr_h0 + min(snr.*snr_SR, snr.*snr_RD);
    
    Pout_HD_DF(isnr) = mean(snr_e2e_HD_DF < snr_th_HD );
    
    EC_HD_DF(isnr) = mean(log2(1 + snr_e2e_HD_DF))/2;
    
    % FD-AF
    snr_e2e_FD_AF =   snr*snr_h0 + (snr.*snr_SR.*snr.*snr_RD) ...
        ./(snr.*snr_SR + (snr.*snr_RD+1).*(snr.*snr_LI+1) );
    
    Pout_FD_AF(isnr) = mean( N_RIS.*snr_e2e_FD_AF < snr_th_FD );
    
    EC_FD_AF(isnr) = mean(log2(1 + N_RIS*snr_e2e_FD_AF));
    
    % FD-DF
    snr_e2e_FD_DF = snr*snr_h0 + min(snr.*snr_SR./(snr.*snr_LI + 1), snr.*snr_RD);
    
    Pout_FD_DF(isnr) = mean(N_RIS.*snr_e2e_FD_DF < snr_th_FD );
    
    EC_FD_DF(isnr) = mean(log2(1 + N_RIS*snr_e2e_FD_DF));
    
end

figure;
semilogy(snrdB,P_out_CCT_Gamma_sim,'rs-'); hold on;
semilogy(snrdB,P_out_OST_Gamma_sim,'bo-'); hold on;
semilogy(snrdB,Pout_nonRIS,'gd-'); hold on;
semilogy(snrdB,Pout_HD_AF,'c+-'); hold on;
semilogy(snrdB,Pout_HD_DF,'c>-'); hold on;
semilogy(snrdB,Pout_FD_AF,'m*-'); hold on;
semilogy(snrdB,Pout_FD_DF,'m^-'); hold on;


xlabel('Transmit SNR (dB)');
ylabel('Outage Probability');
legend('CCT (sim.)',...
    'OST (sim.)',...
    'Non-RIS',...
    'HD-AF',...
    'HD-DF',...
    'FD-AF',...
    'FD-DF',...
    'Location','Best');
axis([-Inf Inf 10^(-5) 10^(0)]);

% x0 = 100; y0 = 100; width = 400; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====

%% Ergodic Capacity

figure;
plot(snrdB,EC_CCT_sim, 'rs-'); hold on;
plot(snrdB,EC_OST_sim, 'bo-'); hold on;
plot(snrdB, EC_nonRIS, 'gd-'); hold on;
plot(snrdB, EC_HD_AF, 'c+-'); hold on;
plot(snrdB, EC_HD_DF, 'c>-'); hold on;
plot(snrdB, EC_FD_AF, 'm*-'); hold on;
plot(snrdB, EC_FD_DF, 'm^-'); hold on;

xlabel('Transmit SNR (dB)');
ylabel('Ergodic Capacity (bits/sec/Hz');
legend('CCT (sim.)',...
    'OST (sim.)',...
    'Non-RIS',...
    'HD-AF',...
    'HD-DF',...
    'FD-AF',...
    'FD-DF',...
    'Location','nw');

% x0 = 100; y0 = 100; width = 400; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====

%% Network Topology

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
