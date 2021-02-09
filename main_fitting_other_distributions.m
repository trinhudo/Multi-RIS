clear all
% close all
%
sim_times = 1e6; % number of simulation trails
%
snrdB = -10:.5:10; % average transmit SNR
r_th = 1; % predefined data rate bits/s/Hz
snr_th = 2^(r_th)-1; % predefined SNR threshold
%
N_RIS = 3; % number of RISs
L = [25 25 25]; % number of elements at each RIS
%
kappa_nl = 1;
% m_h = linspace(1, 3, N_RIS)';
% m_g = linspace(1, 3, N_RIS)';
m0 = 1.5;
m_h = [3.1, 2.7, 2.3]';
m_g = [2.3, 2.7, 3.1]';
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

%%

%-----Monte-Carlo simulation
h_0 = random('Naka',m0,Omg0,[1,sim_times]);

V_n = zeros(N_RIS,sim_times);
for n = 1:N_RIS
    for l = 1:L(n)
        %
        h_nl = random('Naka',m_h(n),Omgh(n),[1,sim_times]);
        g_nl = random('Naka',m_g(n),Omgg(n),[1,sim_times]);
        U_nl = kappa_nl * h_nl .* g_nl;
        %
        V_n(n,:) = V_n(n,:) + U_nl;
    end
end
T_cct = sum(V_n,1);
Z_cct = h_0 + T_cct;
Z2_cct = Z_cct.^2;

V_M_ost = max(V_n,[],1);
R_ost = h_0 + V_M_ost;
R2_ost = R_ost.^2;
%=====DONE

%% Fitting to different distributions

pd_Gamma = fitdist(Z_cct', 'Gamma');
a_Gamma = pd_Gamma.a;
b_Gamma = pd_Gamma.b;
data_Gamma = gamrnd(a_Gamma, b_Gamma, [1, sim_times]);

pd_Lognormal = fitdist(Z_cct', 'Lognormal'); % Z^2 fit to LogNormal
mu_Lognormal = pd_Lognormal.mu;
sigma_Lognormal = pd_Lognormal.sigma;
data_Lognormal = lognrnd(mu_Lognormal, sigma_Lognormal, [1, sim_times]);

pd_Normal = fitdist(Z_cct', 'Normal');
mu_Normal = pd_Normal.mu;
sigma_Normal = pd_Normal.sigma;
data_Normal = normrnd(mu_Normal, sigma_Normal, [1, sim_times]);

pd_Burr = fitdist(Z_cct', 'Burr');
alpha_Burr = pd_Burr.alpha;
c_Burr = pd_Burr.c;
k_Burr = pd_Burr.k;
data_Burr = random('burr', alpha_Burr, c_Burr, k_Burr, 1, sim_times);

pd_Weibull = fitdist(Z_cct', 'Weibull');
A_Weibull = pd_Weibull.A;
B_Weibull = pd_Weibull.B;
data_Weibull = wblrnd(A_Weibull, B_Weibull, [1, sim_times]);

%% Testing in terms of OUTAGE PROBABILITY

%-----
for isnr = 1:length(snrdB)
    snr = 10^(snrdB(isnr)/10);
    %
    Pout_CCT_sim(isnr) = mean(snr*Z2_cct < snr_th);
    Pout_Gamma(isnr) = mean(snr*(data_Gamma.^2) < snr_th);
    Pout_Lognormal(isnr) = mean(snr*(data_Lognormal.^2) < snr_th);
    Pout_Normal(isnr) = mean(snr*(data_Normal.^2) < snr_th);
    Pout_Weilbull(isnr) = mean(snr*(data_Weibull.^2) < snr_th);
    Pout_Burr(isnr) = mean(snr*(data_Burr.^2) < snr_th);
    %
    fprintf('Outage Probability, SNR = %d \n', snrdB(isnr));
end

figure;
semilogy(snrdB,Pout_CCT_sim,'k-', 'linewidth', 2); hold on;
semilogy(snrdB,Pout_Gamma,'r-', 'linewidth', 2); hold on;
semilogy(snrdB,Pout_Lognormal,'g--', 'linewidth', 2); hold on;
semilogy(snrdB,Pout_Normal,'b:', 'linewidth', 2); hold on;
semilogy(snrdB,Pout_Burr,'c-.', 'linewidth', 2); hold on;
semilogy(snrdB,Pout_Weilbull,'m-.', 'linewidth', 2); hold on;


xlabel('Transmit SNR (dB)');
ylabel('Outage Probability');
legend('True (Sim)',...
    'Gamma',...
    'LogNormal',...
    'Normal',...
    'Burr',...
    'Weibull',...
    'Location','Best');
axis([-Inf Inf 10^(-5) 10^(0)]);

figure;
[f_true, x_true] = ksdensity(Z_cct);
plot(x_true, f_true, 'k-', 'linewidth', 2); 
hold on

[f_Gamma, x_Gamma] = ksdensity(data_Gamma);
plot(x_Gamma, f_Gamma, 'r-', 'linewidth', 2); 
hold on

[f_Lognormal, x_Lognormal] = ksdensity(data_Lognormal);
plot(x_Lognormal, f_Lognormal, 'g-', 'linewidth', 2);
hold on

[f_Normal, x_Normal] = ksdensity(data_Normal);
plot(x_Normal, f_Normal, 'b-', 'linewidth', 2);
hold on

[f_Burr, x_Burr] = ksdensity(data_Burr);
plot(x_Burr, f_Burr, 'c-', 'linewidth', 2);
hold on 

[f_Weibull, x_Weibull] = ksdensity(data_Weibull);
plot(x_Weibull, f_Weibull, 'm-', 'linewidth', 2);
hold on 


xlabel('x');
ylabel('Density');
legend('True (Sim)',...
    'Gamma',...
    'LogNormal',...
    'Normal',...
    'Burr',...
    'Weibull',...
    'Location','Best');