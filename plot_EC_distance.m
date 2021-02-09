clear all
close all
%
sim_times = 1e4; % number of simulation trails
%
snrdB = 10; % average transmit SNR
r_th = 1.0; % predefined data rate bits/s/Hz
snr_th = 2^(r_th)-1; % predefined SNR threshold
%
N_RIS = 10; % number of RISs
L = zeros(1, N_RIS); % number of elements at each RIS
%
kappa_nl = 1;
% m_h = linspace(2, 3, N_RIS)';
% m_g = linspace(3, 2, N_RIS)';
m0 = 1.5;
m_h = linspace(2, 3, N_RIS)';
m_g = linspace(3, 2, N_RIS)';
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
x_ris = linspace(x_area_min+x_area_max/N_RIS, x_area_max-x_area_max/N_RIS, N_RIS)'; % [num_ris x 1] vector
y_ris = .7* ones(N_RIS,1)* y_area_max;

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

EC_CCT_sim = zeros(1,N_RIS);

%% Data 


EC_cent_10_75 = [6.45716716384929,5.06838817378597,4.36471880167189,3.98979886809587,3.82496653083392,3.82792768279944,3.99243454557287,4.36627435451065,5.06913849953118,6.45851268924598];
EC_cent_10_150 = [8.41907960409560,6.98974754193470,6.25519546683426,5.85396443796713,5.67601360060527,5.67540190019392,5.85646612051040,6.25445753782963,6.99210796935222,8.41796391734353];

EC_dist_10_75 = 5.56332013630972*ones(1, N_RIS);
EC_dist_10_150 = 7.50319391114484*ones(1, N_RIS);

plot(x_ris, EC_cent_10_75, 'ro-'); hold on
plot(x_ris, EC_dist_10_75, 'b*-'); hold on

plot(x_ris, EC_cent_10_150,'ro-'); hold on
plot(x_ris, EC_dist_10_150, 'b*-'); hold on

xlabel('Distance (m)')
ylabel('Ergodic Capacity (bit/s/Hz)')
legend('Centralized', 'Distributed')
