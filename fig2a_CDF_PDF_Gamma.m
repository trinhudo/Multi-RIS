% Multi-RIS-aided Wireless Systems: Statistical Characterization and Performance Analysis
% Tri Nhu Do, Georges Kaddoum, Thanh Luan Nguyen, Daniel Benevides da Costa, Zygmunt J. Haas
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

%% ERA SCHEME WITH GAMMA DISTRIBUTION

lambda = sqrt(m_h./Omg_h .* m_g./Omg_g) ./ kappa_nl; % lambda_nl

%Working on h0
%The k-th moment of h0
E_h0_k = @(k) gamma(m_0+k/2)/gamma(m_0)*(m_0/Omg_0)^(-k/2);

%CDF of h0
F_h0 = @(x) gammainc(m_0*double(x).^2/Omg_0, m_0, 'lower');

%PDF of h0
f_h0 = @(x) 2*m_0^m_0/gamma(m_0)/Omg_0^m_0*double(x).^(2*m_0-1).*exp(-m_0/Omg_0.*double(x).^2);

%--------------------------------------------------------------------------
%Working on U_nl
%The k-moment of U_nl
E_U_nl_k = @(k,n) lambda(n)^(-k)*gamma(m_h(n)+0.5*k)...
    * gamma(m_g(n)+0.5*k) / gamma(m_h(n)) / gamma(m_g(n));

%Parameter of the approximate Gamma distribution of U_nl
alpha_U= @(n) E_U_nl_k(1,n)^2/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);
beta_U = @(n) E_U_nl_k(1,n)/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);

%PDF of U_nl
f_U_nl = @(x,n) beta_U(n)^alpha_U(n)/gamma(alpha_U(n))...
    * x.^(alpha_U(n)-1) .* exp( -beta_U(n)*x );

%--------------------------------------------------------------------------
%Working on V_n
%The k-moment of V_n
E_V_n_k = @(k,n) gamma(L(n) * alpha_U(n)+k) ...
    / gamma(L(n) * alpha_U(n)) * beta_U(n)^(-k);

%PDF of V_n
f_V_n = @(v,n) vpa(beta_U(n)^(sym(L(n)*alpha_U(n)))/gamma(sym(L(n)*alpha_U(n))))...
    * v.^(L(n)*alpha_U(n)-1) .* exp(-beta_U(n)*v);

%CDF of V_n
F_V_n = @(v,n) gammainc(beta_U(n)*double(v),L(n)*alpha_U(n),'lower');

%--------------------------------------------------------------------------
%Working on T
%The 1st moment of T
E_T_1 = 0;

for nn = 1:N_RIS
    for kk = 1:L(nn)
        E_T_1 = E_T_1 + E_U_nl_k(1,nn);
    end
end

%The 2nd  moment of T
E_T_2 = 0;

for nn = 1:N_RIS
    tmpA = 0;
    for kk = 1:L(nn)
        tmpA = tmpA + E_U_nl_k(1,nn);
    end
    for ii = nn+1:N_RIS
        tmpB = 0;
        for kk = 1:L(ii)
            tmpB = tmpB + E_U_nl_k(1,ii);
        end
        E_T_2 = E_T_2 + 2 * tmpA * tmpB;
    end
end

for nn = 1:N_RIS
    tmpC = 0;
    for kk = 1:L(nn)
        tmpC = tmpC + E_U_nl_k(2,nn);
    end
    tmpD = 0;
    for kk = 1:L(nn)
        for v = (kk+1):L(nn)
            tmpD = tmpD + 2 * E_U_nl_k(1,nn) * E_U_nl_k(1,nn);
        end
    end
    E_T_2 = E_T_2 + tmpC + tmpD;
end

%--------------------------------------------------------------------------
%Working on Z
%The 1st moment of Z
E_Z_1 = E_h0_k(1) + E_T_1;
%The 2nd moment of Z
E_Z_2 = E_h0_k(2) + E_T_2 + 2*E_h0_k(1)*E_T_1;

%Parameter of the approximate Gamma distribution of Z
alpha_Z = E_Z_1^2/(E_Z_2 - E_Z_1^2);
beta_Z = E_Z_1/(E_Z_2 - E_Z_1^2);

%CDF of Z
F_Z_Gamma = @(z) gammainc(z*beta_Z, alpha_Z, 'lower');

%PDF of Z
f_Z_Gamma = @(z) 1/gamma(alpha_Z)*(beta_Z)^alpha_Z...
    * z.^(alpha_Z-1) .* exp( -z*beta_Z );

%CDF of Z^2
F_Z2_Gamma = @(z) F_Z_Gamma(sqrt(z));

%--------------------------------------------------------------------------
%Asymptotic CDF of Z
F_Z_Gamma_asymp = @(z) (z*beta_Z)^alpha_Z/gamma(alpha_Z+1);

%Asymptotic CDF of Z^2
F_Z2_Gamma_asymp= @(z) F_Z_Gamma_asymp(sqrt(z));

%% ORA SCHEME WITH GAMMA DISTRIBUTION

%Working on M_V ( max V_n )
%CDF of V_M
F_M_V = @(x) 1;
for k = 1:N_RIS
    F_M_V = @(x) F_M_V(x) .* F_V_n(x,k);
end

M = 100; %Number of steps in M-staircase approximation

%CDF of R
F_R = @(r) 0;
for m = 1:M
    F_R = @(r) F_R(r)...
        + (F_h0(m/M*r) - F_h0((m-1)/M*r)) .* F_M_V((M-m+1)/M*r);
end

%CDF of R^2
F_R2_Gamma = @(r) F_R(sqrt(r)); %CDF of R^2, e2e SNR of the ORA scheme

%--------------------------------------------------------------------------
% %Asymptotic of R^2
% LAlpha_arr = zeros(1,N_RIS);
% Beta_arr  = zeros(1,N_RIS);
%
% for n = 1:N_RIS
%     LAlpha_arr(n) = L(n)*alpha_U(n);
%     Beta_arr(n) = beta_U(n);
% end
%
% M1 = 1e3; m_arr = 1:M1;
% fM = sum( ((m_arr-1)/M1).^(2*m0-1).*(1-(m_arr-1)/M1).^sum(LAlpha_arr) )/M1;
%
% F_R_Gamma_asymp = @(r) (2*m0)/gamma(m0+1)*(m0/Omg0)^(m0)...
%     * prod( Beta_arr.^(LAlpha_arr) ./ gamma(LAlpha_arr+1) .* r.^LAlpha_arr )...
%     * r.^(2*m0) * fM;
% F_R2_asymp = @(r) F_R_Gamma_asymp(sqrt(r));

%% NETWORK TOPOLOGY

subplot(3,1,1);
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

xlabel('$d_{\rm SD}$ (m)', 'Interpreter', 'Latex')
ylabel('$H$ (m)', 'Interpreter', 'Latex')
axis([x_area_min x_area_max y_area_min y_area_max])
legend('$\rm S$', '$\rm D$', '$\mathrm{R}_i$',...
    'Interpreter', 'Latex',...
    'Location', 'se')

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%% CDFs and PDFs in THE ERA SCHEME

%CDF of Z
% figure(2);
subplot(3,1,2);
[y, x] = ecdf(Z_ERA); hold on;
domain_Z = linspace(0, max(x), 30);

plot(x, y); hold on;
plot(domain_Z, F_Z_Gamma(domain_Z), '.', 'markersize', 10); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('CDF','Interpreter', 'Latex')
legend('True',...
    'Approx.',...
    'location', 'se',...
    'Interpreter', 'Latex');

% x0 = 100; y0 = 100; width = 300; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
set(gca,'fontsize',13);

%--------------------------------------------------------------------------
%Plot PDF of Z

% %non-symbolic PDF of Z and Z^2
% f_Z_Gamma = @(z) 1/gamma(alpha_Z)*(beta_Z)^alpha_Z...
%     * z.^(alpha_Z-1) .* exp( -z*beta_Z );
% f_Z2_Gamma = @(z) 1./(2*sqrt(z)) .* f_Z_Gamma(sqrt(z));

subplot(3,1,3);

number_of_bins = 30;
histogram(Z_ERA, number_of_bins, 'normalization', 'pdf'); hold on;
% histfit(Z_ERA, number_of_bins, 'gamma'); hold on;

%Symbolic PDF of Z
syms sbl_az sbl_bz sbl_z
symbolic_f_Z_Gamma(sbl_az,sbl_bz,sbl_z) = 1/gamma(sbl_az)*(sbl_bz)^sbl_az ...
    * sbl_z^(sbl_az-1) * exp( - sbl_z*sbl_bz );

plot(domain_Z, double(vpa(symbolic_f_Z_Gamma(alpha_Z, beta_Z, domain_Z))),...
    'linewidth', 1); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('PDF', 'Interpreter', 'Latex')
legend('True',...
    'Approx.',...
    'location', 'ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%% CDFs and PDFs in THE ORA SCHEME

%Plot CDF of R
% figure(2);
subplot(3,1,2);
[y, x] = ecdf(R_ORA); hold on;
domain_R = linspace(0, max(x), 30);

plot(x, y); hold on;
plot(domain_R, F_R(domain_R), '.', 'markersize', 10); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('CDF', 'Interpreter', 'Latex')
legend('True',...
    'Approx.',...
    'location', 'se',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);
%--------------------------------------------------------------------------
%Plot PDF of R
%CDF of M_v
f_M_V = @(x) 0;

for nn = 1:N_RIS
    func_tmp = @(x) 1;
    for t = 1:N_RIS
        if (nn ~= t)
            func_tmp = @(x) func_tmp(x) .* F_V_n(x,t);
        end
    end
    f_M_V = @(x) f_M_V(x) + f_V_n(x,nn) .* func_tmp(x);
end

M = 100; %Number of steps in M-staircase approximation

%CDF of R
f_R = @(r) 0;

for m = 1:M
    f_R = @(r) f_R(r)...
        + ((m/M)*f_h0(m/M*r) - ((m-1)/M)*f_h0((m-1)/M*r))...
        .*F_M_V((M-m+1)/M*r)...
        + (F_h0(m/M*r)- F_h0((m-1)/M*r))...
        .*((M-m+1)/M).*f_M_V((M-m+1)/M*r);
end

%CDF of R^2
f_R2 = @(r) 1./(2*sqrt(r)).*f_R(sqrt(r));

% figure(3);
subplot(3,1,3);
histogram(R_ORA, number_of_bins, 'normalization', 'pdf'); hold on;
plot(domain_R, double(vpa(f_R(sym(domain_R)))), 'linewidth', 1.5); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('PDF', 'Interpreter', 'Latex')
legend('True ',...
    'Approx.',...
    'location', 'ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%% OUTAGE PROBABILITY

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); %i.e., 10^(SNRdB/10)
    
    OP_SISO(idx) = mean(avgSNR*SNR_h0 < SNR_th);
    
    %ERA scheme
    OP_ERA_Gamma_sim(idx) = mean(avgSNR*Z2_ERA < SNR_th);
    OP_ERA_Gamma_ana(idx) = F_Z2_Gamma(SNR_th/avgSNR);
    
    %ORA scheme
    OP_ORA_Gamma_sim(idx) = mean(avgSNR*R2_ORA < SNR_th);
    OP_ORA_Gamma_ana(idx) = F_R2_Gamma(SNR_th/avgSNR);
    
    fprintf('Gamma, Outage probability, SNR = %d \n', round(SNRdB(idx)));
end

figure;
semilogy(P_S_dB, OP_SISO, 'gd-'); hold on;

semilogy(P_S_dB, OP_ERA_Gamma_sim, 'r*:'); hold on;
semilogy(P_S_dB, OP_ERA_Gamma_ana, 'r-'); hold on;

semilogy(P_S_dB, OP_ORA_Gamma_sim, 'bo:'); hold on;
semilogy(P_S_dB, OP_ORA_Gamma_ana, 'b-'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Outage probability, $P_{\rm out}$', 'Interpreter', 'Latex');
legend('SISO (sim.)',...
    'ERA (sim.)', ...
    'ERA (ana. with Gamma)',...
    'ORA (sim.)', ...
    'ORA (ana. with Gamma)',...
    'Location','SW',...
    'Interpreter', 'Latex');
axis([-Inf Inf 10^(-5) 10^(0)]);

%% ERGODIC CAPACITY

%EC of ERA scheme with Gamma distribution
%Inputs: az = alpha_Z, bz = beta_Z, rho = snr

syms az bz rho
EC_ERA_Gamma(az, bz, rho) = 1/gamma(az)/log(2)*2^(az-1)/sqrt(pi)...
    * meijerG([0], [1/2,1], [az/2,az/2+1/2,0,0,1/2], [], (bz/2)^2/rho);

%--------------------------------------------------------------------------

func_EC_ORA_Gamma = @(x,c) (1/log(2)).*(1./(1+x).*(1 - F_R(sqrt(x./c))));

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); %10^(SNRdB(idx)/10)
    
    EC_SISO(idx) = mean(log2(1 + avgSNR*SNR_h0));
    
    EC_ERA_Gamma_sim(idx) = mean(log2(1+avgSNR*Z2_ERA));
    EC_ERA_Gamma_ana(idx) = double(vpa(EC_ERA_Gamma(alpha_Z, beta_Z, avgSNR)));
    
    EC_ORA_Gamma_sim(idx) = mean(log2(1+avgSNR*R2_ORA));
    EC_ORA_Gamma_ana(idx) = integral(@(x) func_EC_ORA_Gamma(x,avgSNR), 0, Inf);
    
    fprintf('Gamma, Ergodic capacity, SNR = %d \n', round(SNRdB(idx)));
end

figure;
plot(P_S_dB, EC_SISO, 'go-'); hold on;
plot(P_S_dB, EC_ERA_Gamma_sim, 'ro'); hold on;
plot(P_S_dB, EC_ERA_Gamma_ana, 'r-'); hold on;

plot(P_S_dB, EC_ORA_Gamma_sim, 'bs:'); hold on;
plot(P_S_dB, EC_ORA_Gamma_ana, 'b-'); hold on;

xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [b/s/Hz]', 'Interpreter', 'Latex');
legend('SISO (sim.)',...
    'ERA (sim.)',...
    'ERA (ana. with Gamma)',...
    'ORA (sim.)',...
    'ORA (analytical)',...
    'Interpreter', 'Latex',...
    'Location','NW');

%% SAVE DATA

save('data_Gamma_settingL1.mat',...
    'OP_SISO',...
    'OP_ERA_Gamma_sim',...
    'OP_ERA_Gamma_ana',...
    'OP_ORA_Gamma_sim',...
    'OP_ORA_Gamma_ana',...
    'EC_SISO',...
    'EC_ERA_Gamma_sim',...
    'EC_ERA_Gamma_ana',...
    'EC_ORA_Gamma_sim',...
    'EC_ORA_Gamma_ana')
