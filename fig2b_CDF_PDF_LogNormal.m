% Multi-RIS-aided Wireless Systems: Statistical Characterization and Performance Analysis
% Tri Nhu Do, Georges Kaddoum, Thanh Luan Nguyen, Daniel Benevides da Costa, Zygmunt J. Haas
% https://arxiv.org/abs/2104.01912
% Version: 2021-04-05

%% Simulation parameters

clear;
close all;

sim_trials = 1e5; %Number of simulation trails

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

% %Loacation setting D1
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

%% ANALYSIS OF THE ERA SCHEME

lambda = sqrt(m_h./Omg_h .* m_g./Omg_g) ./ kappa_nl; % lambda_nl

% Working on h0
% The k-th moment of h0
E_h0_k = @(k) gamma(m_0+k/2)/gamma(m_0)*(m_0/Omg_0)^(-k/2);

F_h0 = @(x) gammainc(m_0*x.^2/Omg_0, m_0, 'lower');

%--------------------------------------------------------------------------
% Working on U_nl
% The k-moment of U_nl
E_U_nl_k = @(k, n) lambda(n)^(-k)*gamma(m_h(n)+0.5*k)...
    * gamma(m_g(n)+0.5*k) / gamma(m_h(n)) / gamma(m_g(n));

% Parameter of the approximate Gamma distribution of U_nl
alpha_U = @(n) E_U_nl_k(1, n)^2/(E_U_nl_k(2, n)-E_U_nl_k(1, n)^2);
beta_U = @(n) E_U_nl_k(1, n)/(E_U_nl_k(2, n)-E_U_nl_k(1, n)^2);

%--------------------------------------------------------------------------
% Working on V_n
% The k-moment of V_n
E_V_n_k = @(k, n) gamma(L(n) * alpha_U(n)+k) ...
    / gamma(L(n) * alpha_U(n)) * beta_U(n)^(-k);

%--------------------------------------------------------------------------
%Working on T
%The 1st moment of T
E_T1 = 0;
for n = 1:N_RIS
    for l = 1:L(n)
        E_T1 = E_T1 + E_U_nl_k(1, n);
    end
end

%The 2nd  moment of T
E_T2 = 0;
for n = 1:N_RIS
    tmpA = 0;
    for l = 1:L(n)
        tmpA = tmpA + E_U_nl_k(1, n);
    end
    for ii = n+1:N_RIS
        tmpB = 0;
        for l = 1:L(ii)
            tmpB = tmpB + E_U_nl_k(1, ii);
        end
        E_T2 = E_T2 + 2 * tmpA * tmpB;
    end
end

for n = 1:N_RIS
    tmpC = 0;
    for l = 1:L(n)
        tmpC = tmpC + E_U_nl_k(2, n);
    end
    tmpD = 0;
    for l = 1:L(n)
        for v = (l+1):L(n)
            tmpD = tmpD + 2 * E_U_nl_k(1, n) * E_U_nl_k(1, n);
        end
    end
    E_T2 = E_T2 + tmpC + tmpD;
end

E_Z = E_h0_k(1) + E_T1;
E_Z1 = E_h0_k(1) + E_T1;
E_Z2 = E_h0_k(2) + E_T2 + 2 * E_h0_k(1) * E_T1;

%--------------------------------------------------------------------------
% Fit Z2_ERA to Log-Normal
E_Z4 = 0; % the 2nd moment of Z^2, i.e., E[(Z^2)^2]
for k = 0:4
    %
    E_Tk = 0; % E[T^k] % the k-th momemt of T >> CAN BE USE IN ERA ???
    [cases_T, indT_mat] = nsumk(N_RIS, k);
    for icaseT = 1:cases_T
        indT_arr = indT_mat(icaseT, :);
        %
        tmpT = 1;
        for t = 1:N_RIS
            tmpT = tmpT * E_V_n_k(indT_arr(t), t);
        end
        %
        E_Tk = E_Tk + factorial(k)/prod( factorial(indT_arr) ) * tmpT;
        %
    end
    %
    E_Z4 = E_Z4 + nchoosek(4, k)*E_h0_k(4-k)*E_Tk;
end

nu_ERA_LN = log( E_Z2^2/sqrt(E_Z4) ); % nu_Z
zeta2_ERA_LN = log( E_Z4/E_Z2^2 ); % zeta_Z^2

nu_ERA_new = log(E_Z^2 / sqrt(E_Z2) );
zeta2_ERA_new = log( E_Z2 / (E_Z^2) );

%CDF of Z 
F_Z_ERA_new = ...
    @(x) 1/2 + 1/2*erf( (log(x)-nu_ERA_new)/sqrt(2*zeta2_ERA_new) ); 

%CDF of Z^2
F_Z2_ERA_LN = ...
    @(x) 1/2 + 1/2*erf( (log(x)-nu_ERA_LN)/sqrt(2*zeta2_ERA_LN) ); 


%% ANALYSIS OF THE ORA SCHEME

alpha_U_arr = zeros(1, N_RIS);
beta_U_arr = zeros(1, N_RIS);

for n = 1:N_RIS
    alpha_U_arr(n) = alpha_U(n);
    beta_U_arr(n) = beta_U(n);
end
%
chi_t= @(t) beta_U_arr(t) ./ sum(beta_U_arr);

%--------------------------------------------------------------------------
% Approxiate result, using self-built F_A() function

alpha_U_arr= zeros(1, N_RIS);
beta_U_arr = zeros(1, N_RIS);
for n = 1:N_RIS
    alpha_U_arr(n) = alpha_U(n);
    beta_U_arr(n) = beta_U(n);
end

f_V_n = @(v, n) beta_U(n)^(L(n)*alpha_U(n))/gamma(L(n)*alpha_U(n))...
    * v.^(L(n)*alpha_U(n)-1) .* exp( -beta_U(n)*v );
F_V_n = @(v, n) gammainc(beta_U(n)*v, L(n)*alpha_U(n), 'lower');

f_M_V = @(x) 0;

for n = 1:N_RIS
    func_tmp = @(x) 1;
    for t = 1:N_RIS
        if (n ~= t)
            func_tmp = @(x) func_tmp(x) .* F_V_n(x, t);
        end
    end
    f_M_V = @(x) f_M_V(x) + f_V_n(x, n) .* func_tmp(x);
end

% mu_M_V = zeros(1, 4); % the k-th moment of M_V (k = 1, 2, 3, 4)
% for k = 1:4
%     mu_M_V(k) = integral(@(x) x.^k .* f_M_V(x), 0, 250);
% end

X = sym('X', [1, N_RIS]);
mu_M_V = sym(zeros(1, 4)); % the k-th moment of M_V (k = 1, 2, 3, 4)
for k = 1:4
    for n = 1:N_RIS
        %
        Sn = setdiff(1:N_RIS, n);
        %
        tmp = Lauricella_FA(sum(L.*alpha_U_arr)+k, ones(1, N_RIS-1), L(Sn).*alpha_U_arr(Sn)+1, chi_t(Sn));
        %
        if (tmp > 0)
            mu_M_V(k) = mu_M_V(k) + gamma( sum(X)+k ) / gamma(X(n))...
                / prod( gamma(X(Sn)+1) ) * sym(tmp);
        end
    end
    mu_M_V(k) = vpa(subs(sum(beta_U_arr)^(-k) * prod( chi_t(1:N_RIS).^(X) ) * mu_M_V(k), X, L.*alpha_U_arr));
end
mu_M_V = double(mu_M_V);


E_R = 0; % E[R]
for k = 0:1
    if k >= 1
        E_R = E_R + nchoosek(1, k) * E_h0_k(1-k) * mu_M_V(k);
    else
        E_R = E_R + E_h0_k(1);
    end
end

E_R_2 = 0; % E[R^2] by using R^2 expressions, not R
for k = 0:2
    if k >= 1
        E_R_2 = E_R_2 + nchoosek(2, k) * E_h0_k(2-k) * mu_M_V(k);
    else
        E_R_2 = E_R_2 + E_h0_k(2);
    end
end

E_R_4 = 0; % E[(R^2)^2] by using R^2 expressions, not R
for k = 0:4
    if k >= 1
        E_R_4 = E_R_4 + nchoosek(4, k) * E_h0_k(4-k) * mu_M_V(k);
    else
        E_R_4 = E_R_4 + E_h0_k(4);
    end
end

nu_R_ORA_LN = log( E_R^2/sqrt(E_R_2) );
zeta2_R_ORA_LN = log( E_R_2/E_R^2 ); % zeta^2 of R

nu_R2_ORA_LN = log( E_R_2^2/sqrt(E_R_4) );
zeta2_R2_ORA_LN = log( E_R_4/E_R_2^2 ); % zeta^2 of R^2

F_R_ORA_LN = @(x) 1/2 + ...
    1/2*erf( (log(x)-nu_R_ORA_LN)/sqrt(2*zeta2_R_ORA_LN) ); % CDF of R

F_R2_ORA_LN = @(x) 1/2 + ...
    1/2*erf( (log(x)-nu_R2_ORA_LN)/sqrt(2*zeta2_R2_ORA_LN) ); % CDF of R^2

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

%% CDFs and PDF in the ERA scheme

%CDF of Z based on the LogNormal distribution

subplot(3, 1, 2);

[y, x] = ecdf(Z_ERA); hold on;
domain_Z = linspace(0, max(x), 30);

plot(x, y); hold on;
plot(domain_Z, F_Z_ERA_new(domain_Z), '.', 'markersize', 10); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('CDF','Interpreter', 'Latex')
legend('True',...
    'Approx.',...
    'location', 'se',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
set(gca,'fontsize',13);

% -------------------------------------------------------------------------
%PDF of Z based on the  LogNormal distribution

subplot(3, 1, 3);

f_Z_ERA_new = @(x) 1./(x .* sqrt(2*pi*zeta2_ERA_new) )...
    .* exp( -(log(x) - nu_ERA_new).^2./(2*zeta2_ERA_new) ); % PDF of Z

f_Z2_ERA_LN = @(x) 1./(x .* sqrt(2*pi*zeta2_ERA_LN) )...
    .* exp( -(log(x) - nu_ERA_LN).^2./(2*zeta2_ERA_LN) ); % PDF of Z^2

number_of_bins = 30;
% histogram(Z2_ERA, number_of_bins, 'normalization', 'pdf'); hold on;
histogram(Z_ERA, number_of_bins, 'normalization', 'pdf'); hold on;

plot(domain_Z, double(vpa(f_Z_ERA_new(sym(domain_Z)))), 'linewidth', 1); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('PDF', 'Interpreter', 'Latex')

legend('True',...
    'Approx.',...
    'location', 'ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%% CDFs and PDF in the ORA scheme

%CDF of R based on the LogNormal Distribution
subplot(3, 1, 2);

[y, x] = ecdf(R_ORA); hold on;
domain_R = linspace(0, max(x), 30);
plot(x, y); hold on;
plot(domain_R, F_R_ORA_LN(domain_R), '.', 'markersize', 10); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('CDF', 'Interpreter', 'Latex')

legend('True',...
    'Approx.',...
    'location', 'se',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

%--------------------------------------------------------------------------
%PDF of R based on the LogNormal Distribution

subplot(3, 1, 3);

f_R_ORA_LN = @(x) 1./(x .* sqrt(2*pi*zeta2_R_ORA_LN) )...
    .* exp( -(log(x) - nu_R_ORA_LN).^2./(2*zeta2_R_ORA_LN) ); % CDF of R

%PDF of R^2
f_R2_ORA_LN = @(x) 1./(x .* sqrt(2*pi*zeta2_ORA_LN) )...
    .* exp( -(log(x) - nu_ORA_LN).^2./(2*zeta2_ORA_LN) );

number_of_bins = 30;

histogram(R_ORA, number_of_bins, 'normalization', 'pdf'); hold on;
plot(domain_R, double(vpa(f_R_ORA_LN(sym(domain_R)))), 'linewidth', 1.5); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('PDF', 'Interpreter', 'Latex')
legend('True ',...
    'Approx.',...
    'location', 'ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);

% %% OUTAGE PROBABILITY
% 
% for iSNR = 1:length(SNRdB)
%     SNR = 10^(SNRdB(iSNR)/10);
%     %
%     OP_ERA_LN_sim(iSNR) = mean(SNR*Z2_ERA < SNR_th);
%     OP_ERA_LN_ana(iSNR) = F_Z_ERA_new(sqrt(SNR_th/SNR)); % F_Z (sqrt(x))
%     %
%     OP_ORA_LN_sim(iSNR) = mean(SNR*R2_ORA < SNR_th);
%     OP_ORA_LN_ana(iSNR) = F_R_ORA_LN(sqrt(SNR_th/SNR));
%     fprintf('LogNormal, Outage Probability, SNR = %d \n', SNRdB(iSNR));
% end
% 
% figure;
% semilogy(P_S_dB, OP_ERA_LN_sim, 'r+:'); hold on;
% semilogy(P_S_dB, OP_ERA_LN_ana, 'r-'); hold on;
% semilogy(P_S_dB, OP_ORA_LN_sim, 'bo:'); hold on;
% semilogy(P_S_dB, OP_ORA_LN_ana, 'b-'); hold on;
% 
% xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
% ylabel('Outage Probability, $P_{\rm out}$', 'Interpreter', 'Latex');
% legend('ERA (sim.)', ...
%     'ERA (ana. with LogNormal)',...
%     'ORA (sim.)', ...
%     'ORA (ana. with LogNormal)',...
%     'Location','SW',...
%     'Interpreter', 'Latex');
% axis([-Inf Inf 10^(-5) 10^(0)]);
% 
% %% ERGODIC CAPACITY
% 
% a1 =  0.9999964239;
% a2 = -0.4998741238;
% a3 =  0.3317990258;
% a4 = -0.2407338084;
% a5 =  0.1676540711;
% a6 = -0.0953293897;
% a7 =  0.0360884937;
% a8 = -0.0064535442;
% a_arr = [a1, a2, a3, a4, a5, a6, a7, a8];
% 
% syms aXi bXi Xi(aXi, bXi)
% Xi(aXi, bXi) = 0;
% 
% for k = 1:8
%     Xi(aXi, bXi) = Xi(aXi, bXi) + exp(-bXi^2)/2 * a_arr(k)...
%         * exp((k/(2*aXi) + bXi)^2) * erfc(k/(2*aXi) + bXi);
% end
% 
% syms zeta2 nu
% EC_LN(zeta2, nu) = ...
%     Xi(1/sqrt(2*zeta2),  nu/sqrt(2*zeta2))...
%     + Xi(1/sqrt(2*zeta2), -nu/sqrt(2*zeta2))...
%     + sqrt(zeta2/2/pi) * exp(-nu^2/(2*zeta2))...
%     + nu/2*erfc(-nu/sqrt(2*zeta2));
% 
% for iSNR = 1:length(SNRdB)
%     SNR = 10^(SNRdB(iSNR)/10);
%     %
%     EC_ERA_LN_sim(iSNR) = mean(log2(1 + SNR*Z2_ERA));
%     EC_ERA_LN_ana(iSNR) = double(vpa(EC_LN(zeta2_ERA_LN, log(SNR) + nu_ERA_LN)))/log(2);
%     
%     EC_ORA_LN_sim(iSNR) = mean(log2(1+SNR*R2_ORA));
%     EC_ORA_LN_ana(iSNR) = double(vpa(EC_LN(zeta2_R2_ORA_LN, log(SNR) + nu_R2_ORA_LN)))/log(2);
%     
%     fprintf('LogNormal, Ergodic Capacity , SNR = %d \n', round(SNRdB(iSNR)));
% end
% 
% figure;
% plot(P_S_dB, EC_ERA_LN_sim, 'r+:'); hold on;
% plot(P_S_dB, EC_ERA_LN_ana, 'r-'); hold on;
% plot(P_S_dB, EC_ORA_LN_sim, 'bo:'); hold on;
% plot(P_S_dB, EC_ORA_LN_ana, 'b-'); hold on;
% 
% xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
% ylabel('Ergodic Capacity [b/s/Hz]', 'Interpreter', 'Latex');
% legend('ERA (sim.)','ERA (ana. with LogNormal)',...
%     'ORA (sim.)', ...
%     'ORA (ana. with LogNormal)',...
%     'Interpreter', 'Latex',...
%     'Location','NW');


%% SAVE DATA

% save('data_LogNormal_settingL1.mat', ...
%     'OP_ERA_LN_sim', ...
%     'OP_ERA_LN_ana', ...
%     'OP_ORA_LN_sim', ...
%     'OP_ORA_LN_ana', ...
%     'EC_ERA_LN_sim', ...
%     'EC_ERA_LN_ana', ...
%     'EC_ORA_LN_sim', ...
%     'EC_ORA_LN_ana')
