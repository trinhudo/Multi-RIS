% This script is for single-RIS-aided wireless system

%% Simulation 

clear
close all

sim_times = 1e5; % Number of simulation trails

R_th = 1; % Predefined target spectral efficiency [b/s/Hz]
SNR_th = 2^R_th-1; % Predefined SNR threshold

L = 64; % Number of elements at each RIS

kappa_nl = 1; % Amplitude reflection coefficient

% --------------------------------------------------------------------------

% Network area
x_area_min = 0;
x_area_max = 100; % in meters
y_area_min = 0;
y_area_max = 10;

% Source location
x_source = x_area_min;
y_source = y_area_min;

% Destination location
x_des = x_area_max;
y_des = y_area_min;

% Random location setting
x_RIS = x_area_min + (x_area_max-x_area_min)*rand(1, 1); % [num_RIS x 1] vector
y_RIS = y_area_min + (y_area_max-y_area_min)*rand(1, 1);

% Compute location of nodes
pos_source = [x_source, y_source];
pos_des = [x_des, y_des];
pos_RIS = [x_RIS, y_RIS]; % [num_RIS x 2] matrix

% Compute distances
d_sr = sqrt(sum((pos_source - pos_RIS).^2 , 2)); % [num_RIS x 1] vector
d_rd = sqrt(sum((pos_RIS - pos_des).^2 , 2));
d_sd = sqrt(sum((pos_source - pos_des).^2 , 2));

% --------------------------------------------------------------------------
% Path-loss model
% Carrier frequency (in GHz)
fc = 3; % GHz

% 3GPP Urban Micro in 3GPP TS 36.814, Mar. 2010.
% Note that x is measured in meter

% NLoS path-loss component based on distance
pathloss_NLOS = @(x) db2pow(-22.7 - 26*log10(fc) - 36.7*log10(x));

antenna_gain_S = db2pow(5); % Source antenna gain, dBi
antenna_gain_RIS = db2pow(5); % Gain of each element of a RIS, dBi
antenna_gain_D = db2pow(0); % Destination antenna gain, dBi

% --------------------------------------------------------------------------
% Noise power and Transmit power P_S
% Bandwidth
BW = 10e6; % 10 MHz

% Noise figure (in dB)
noiseFiguredB = 10;

% Compute the noise power in dBm
sigma2dBm = -174 + 10*log10(BW) + noiseFiguredB; % -94 dBm
sigma2 = db2pow(sigma2dBm);

P_S_dB = -5:2:25; % Transmit power of the source, dBm, e.g., 200mW = 23dBm

SNRdB = P_S_dB - sigma2dBm; % Average transmit SNR, dB = dBm - dBm

%% SIMULATION

% Nakagami scale parameter
m_0 = 2.5 + randn; % S->D, scale parameter, heuristic setting
m_h = 2.5 + randn; % S->R
m_g = 2.5 + randn; % R->D

% Nakagami spread parameter
Omega_0 = 1; % Normalized spread parameter of S->D link
Omega_h = 1; % Normalized spread parameter of S->RIS link
Omega_g = 1; % Normalized spread parameter of RIS->D link

% Path-loss
path_loss_0 = pathloss_NLOS(d_sd)*antenna_gain_S; % S->D link

path_loss_h = pathloss_NLOS(d_sr) * ...
    antenna_gain_S*antenna_gain_RIS*L; % Source -> RIS

path_loss_g = pathloss_NLOS(d_rd) * ...
    antenna_gain_RIS*L*antenna_gain_D; % RIS -> Des

% phase of channels
phase_h_SD = 2*pi*rand(1, sim_times); % domain [0,2pi)
phase_h_SR = 2*pi*rand(L, sim_times); % domain [0,2pi)
phase_g_RD = 2*pi*rand(L, sim_times); % domain [0,2pi)

% Channel modeling
h_SD = sqrt(path_loss_0) * ...
    random('Naka', m_0, Omega_0, [1, sim_times]) .* ...
    exp(1i*phase_h_SD);

h_SR = sqrt(path_loss_h) .* ...
    random('Naka', m_h, Omega_h, [L, sim_times]) .* ...
    exp(1i*phase_h_SR);

g_RD = sqrt(path_loss_g) .* ...
    random('Naka', m_g, Omega_g, [L, sim_times]) .* ...
    exp(1i*phase_g_RD);

% Phase-shift configuration

for ss = 1:sim_times % loop over simulation trials
    for ll = 1:L % loop over each elements of the RIS
        % unknown domain phase-shift
        phase_shift_element_temp(ll) = ...
            phase_h_SD(ss) - phase_h_SR(ll,ss) - phase_g_RD(ll,ss); 
        % convert to domain of [0, 2pi)
        phase_shift_element(ll) = wrapTo2Pi(phase_shift_element_temp(ll)); 
        phase_shift_vector(ll) = exp(1i*phase_shift_element(ll));
    end
    phase_shift_matrix = kappa_nl .* diag(phase_shift_vector);
    
    % e2e channel coefficient (complex number, not magnitude)
    h_e2e(:,ss) = h_SD(:,ss) + ...
        h_SR(:,ss).' * phase_shift_matrix * g_RD(:,ss);
end

h_e2e_magnitude = abs(h_e2e); % Magnitude of the e2e channel
h_e2e_squared  = abs(h_e2e_magnitude).^2; % Squared magnitude of the e2e channel


%% Analysis based on moment-matching to Gamma distribution

Omg_0 = Omega_0*path_loss_0;
Omg_h = Omega_h*path_loss_h;
Omg_g = Omega_g*path_loss_g;

lambda = sqrt(m_h./Omg_h .* m_g./Omg_g) ./ kappa_nl; % lambda_nl

% Working on h0
% The k-th moment of h0
E_h0_k = @(k) gamma(m_0+k/2)/gamma(m_0)*(m_0/Omg_0)^(-k/2);

% CDF of h0
F_h0 = @(x) gammainc(m_0*double(x).^2/Omg_0, m_0, 'lower');

% PDF of h0
f_h0 = @(x) 2*m_0^m_0/gamma(m_0)/Omg_0^m_0*double(x).^(2*m_0-1) .* ...
    exp(-m_0/Omg_0.*double(x).^2);

% -------------------------------------------------------------------------
% Working on U_nl
% The k-moment of U_nl
E_U_nl_k = @(k,n) lambda(n)^(-k)*gamma(m_h(n)+0.5*k)...
    * gamma(m_g(n)+0.5*k) / gamma(m_h(n)) / gamma(m_g(n));

% Parameter of the approximate Gamma distribution of U_nl
alpha_U= @(n) E_U_nl_k(1,n)^2/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);
beta_U = @(n) E_U_nl_k(1,n)/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);

% PDF of U_nl
f_U_nl = @(x,n) beta_U(n)^alpha_U(n)/gamma(alpha_U(n))...
    * x.^(alpha_U(n)-1) .* exp( -beta_U(n)*x );

% -------------------------------------------------------------------------
% Working on V_n
% The k-moment of V_n
E_V_n_k = @(k,n) gamma(L(n) * alpha_U(n)+k) ...
    / gamma(L(n) * alpha_U(n)) * beta_U(n)^(-k);

% PDF of V_n
f_V_n = @(v,n) vpa(beta_U(n)^(sym(L(n)*alpha_U(n))) ...
    / gamma(sym(L(n)*alpha_U(n))))...
    * v.^(L(n)*alpha_U(n)-1) .* exp(-beta_U(n)*v);

% CDF of V_n
F_V_n = @(v,n) gammainc(beta_U(n)*double(v),L(n)*alpha_U(n),'lower');

% --------------------------------------------------------------------------
% Working on T
% The 1st moment of T
E_T_1 = 0;

for nn = 1:1
    for kk = 1:L
        E_T_1 = E_T_1 + E_U_nl_k(1,nn);
    end
end

% The 2nd  moment of T
E_T_2 = 0;

for nn = 1:1
    tmpA = 0;
    for kk = 1:L
        tmpA = tmpA + E_U_nl_k(1,nn);
    end
    for ii = nn+1:1
        tmpB = 0;
        for kk = 1:L(ii)
            tmpB = tmpB + E_U_nl_k(1,ii);
        end
        E_T_2 = E_T_2 + 2 * tmpA * tmpB;
    end
end

for nn = 1:1
    tmpC = 0;
    for kk = 1:L
        tmpC = tmpC + E_U_nl_k(2,nn);
    end
    tmpD = 0;
    for kk = 1:L
        for v = (kk+1):L
            tmpD = tmpD + 2 * E_U_nl_k(1,nn) * E_U_nl_k(1,nn);
        end
    end
    E_T_2 = E_T_2 + tmpC + tmpD;
end

% --------------------------------------------------------------------------
% Working on Z
% The 1st moment of Z
E_Z_1 = E_h0_k(1) + E_T_1;
% The 2nd moment of Z
E_Z_2 = E_h0_k(2) + E_T_2 + 2*E_h0_k(1)*E_T_1;

% Parameter of the approximate Gamma distribution of Z
alpha_Z = E_Z_1^2/(E_Z_2 - E_Z_1^2);
beta_Z = E_Z_1/(E_Z_2 - E_Z_1^2);

% CDF of Z
F_Z = @(z) gammainc(z*beta_Z, alpha_Z, 'lower');

% PDF of Z
f_Z = @(z) 1/gamma(alpha_Z)*(beta_Z)^alpha_Z...
    * z.^(alpha_Z-1) .* exp( -z*beta_Z );

% CDF of Z^2
F_Z2 = @(z) F_Z(sqrt(z));

% --------------------------------------------------------------------------
% Asymptotic CDF of Z
F_Z_asymp = @(z) (z*beta_Z)^alpha_Z/gamma(alpha_Z+1);

% Asymptotic CDF of Z^2
F_Z2_asymp= @(z) F_Z_asymp(sqrt(z));


%% NETWORK TOPOLOGY

figure;
scatter(x_source, y_source, 100, 'bo', 'filled')
hold on
scatter(x_des, y_des, 100, 'bd', 'filled')
hold on
scatter(x_RIS, y_RIS, 100, 'rs', 'filled')
hold on

for kk = 1:1
    text(x_RIS(kk)+3, y_RIS(kk), num2str(kk));
    hold on
end

xlabel('$d_{\rm SD}$ (m)', 'Interpreter', 'Latex')
ylabel('$H$ (m)', 'Interpreter', 'Latex')
axis([x_area_min x_area_max y_area_min y_area_max])
legend('$\rm S$', '$\rm D$', '$\mathrm{R}_i$',...
    'Interpreter', 'Latex',...
    'Location', 'se')

set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
set(gca,'fontsize',13);
hold off

%% CDF of Z

figure;

% CDF of Z
subplot(2,1,1);
[y, x] = ecdf(h_e2e_magnitude); hold on;
domain_Z = linspace(0, max(x), 30);

plot(x, y); hold on;
plot(domain_Z, F_Z(domain_Z), '.', 'markersize', 10); hold on;

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

%% PDF of Z

% % non-symbolic PDF of Z and Z^2
% f_Z = @(z) 1/gamma(alpha_Z)*(beta_Z)^alpha_Z...
%     * z.^(alpha_Z-1) .* exp( -z*beta_Z );
% f_Z2 = @(z) 1./(2*sqrt(z)) .* f_Z(sqrt(z));

subplot(2,1,2);

number_of_bins = 30;
histogram(h_e2e_magnitude, number_of_bins, 'normalization', 'pdf'); hold on;
% histfit(h_e2e_magnitude, number_of_bins, 'gamma'); hold on;

% Symbolic PDF of Z
syms sbl_az sbl_bz sbl_z
symbolic_f_Z(sbl_az,sbl_bz,sbl_z) = 1/gamma(sbl_az)*(sbl_bz)^sbl_az ...
    * sbl_z^(sbl_az-1) * exp( - sbl_z*sbl_bz );

plot(domain_Z, double(vpa(symbolic_f_Z(alpha_Z, beta_Z, domain_Z))),...
    'linewidth', 1); hold on;

xlabel('$x$', 'Interpreter', 'Latex')
ylabel('PDF', 'Interpreter', 'Latex')
legend('True',...
    'Approx.',...
    'location', 'ne',...
    'Interpreter', 'Latex');

set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
set(gca,'fontsize',13);

%% OUTAGE PROBABILITY

SNR_h0 = abs(h_SD).^2;

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); % i.e., 10^(SNRdB/10)
    
    OP_SISO(idx) = mean(avgSNR*SNR_h0 < SNR_th);
    
    % RIS scheme
    OP_RIS_sim(idx) = mean(avgSNR*h_e2e_squared < SNR_th);
    OP_RIS_ana(idx) = F_Z2(SNR_th/avgSNR);
    
    fprintf('Outage probability, SNR = % d \n', round(SNRdB(idx)));
end

figure;
semilogy(P_S_dB, OP_SISO, 'b*-'); hold on;

semilogy(P_S_dB, OP_RIS_sim, 'ro:'); hold on;
semilogy(P_S_dB, OP_RIS_ana, 'r-'); hold on;


xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Outage probability, $P_{\rm out}$', 'Interpreter', 'Latex');
legend('Non-RIS (sim.)',...
    'RIS (sim.)', ...
    'RIS (ana.)',...
    'Location','SW',...
    'Interpreter', 'Latex');
axis([-Inf Inf 10^(-5) 10^(0)]);

%% ERGODIC CAPACITY

% EC of RIS scheme with Gamma distribution
% Inputs: az = alpha_Z, bz = beta_Z, rho = snr

syms az bz rho
EC_RIS(az, bz, rho) = 1/gamma(az)/log(2)*2^(az-1)/sqrt(pi)...
    * meijerG(0, [1/2,1], [az/2,az/2+1/2,0,0,1/2], [], (bz/2)^2/rho);

% --------------------------------------------------------------------------

func_EC_ORA = @(x,c) (1/log(2)).*(1./(1+x).*(1 - F_R(sqrt(x./c))));

for idx = 1:length(SNRdB)
    avgSNR = db2pow(SNRdB(idx)); % 10^(SNRdB(idx)/10)
    
    EC_non_RIS(idx) = mean(log2(1 + avgSNR*SNR_h0));
    
    EC_RIS_sim(idx) = mean(log2(1+avgSNR*h_e2e_squared));
    EC_RIS_ana(idx) = double(vpa(EC_RIS(alpha_Z, beta_Z, avgSNR)));
    
    fprintf('Ergodic capacity, SNR = % d \n', round(SNRdB(idx)));
end

figure;
plot(P_S_dB, EC_non_RIS, 'b*-'); hold on;
plot(P_S_dB, EC_RIS_sim, 'ro'); hold on;
plot(P_S_dB, EC_RIS_ana, 'r-'); hold on;


xlabel('$P_{\rm S}$ [dBm]', 'Interpreter', 'Latex');
ylabel('Ergodic capacity [b/s/Hz]', 'Interpreter', 'Latex');
legend('Non-RIS (sim.)',...
    'RIS (sim.)',...
    'RIS (ana.)',...
    'Interpreter', 'Latex',...
    'Location','NW');
