clear all
close all
%
sim_trials = 1e6; % number of simulation trails
%
snrdB = -10:1:15; % average transmit SNR
r_th = 1.0; % predefined data rate bits/s/Hz
snr_th = 2^(r_th)-1; % predefined SNR threshold
%
N_RIS = 3; % number of RISs
L = [30 25 20]; % number of elements at each RIS
%
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
h_0 = random('Naka',m0,Omg0,[1,sim_trials]);

V_n = zeros(N_RIS, sim_trials);
for n = 1:N_RIS
    for l = 1:L(n)
        h_nl = random('Naka', m_h(n), Omgh(n), [1,sim_trials]);
        g_nl = random('Naka', m_g(n), Omgg(n), [1,sim_trials]);
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

%% Analysis of CCT using Gamma distribution

lambda = sqrt( m_g./Omgg.*m_h./Omgh )/kappa_nl; % lambda_nl

% Working on h0
% The k-th moment of h0
E_h0_k = @(k) gamma(m0+k/2)/gamma(m0)*(m0/Omg0)^(-k/2);

F_h0 = @(x) gammainc(m0*double(x).^2/Omg0,m0,'lower'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CDF of h0
f_h0 = @(x) 2*m0^m0/gamma(m0)/Omg0^m0*double(x).^(2*m0-1).*exp(-m0/Omg0.*double(x).^2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % PDF of h0

% Working on U_nl
% The k-moment of U_nl
E_U_nl_k = @(k,n) lambda(n)^(-k)*gamma(m_h(n)+0.5*k)...
    * gamma(m_g(n)+0.5*k) / gamma(m_h(n)) / gamma(m_g(n));

% Parameter of the approximate Gamma distribution of U_nl
alpha_U= @(n) E_U_nl_k(1,n)^2/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);
beta_U = @(n) E_U_nl_k(1,n)/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);

f_U_nl = @(x,n) beta_U(n)^alpha_U(n)/gamma(alpha_U(n))...
    * x.^(alpha_U(n)-1) .* exp( -beta_U(n)*x ); % CDF of U_nl

% Working on V_n
% The k-moment of V_n
E_V_n_k = @(k,n) gamma(L(n) * alpha_U(n)+k) ...
    / gamma(L(n) * alpha_U(n)) * beta_U(n)^(-k);

% Approximate distribution of V_n
f_V_n = @(v,n) vpa(beta_U(n)^(sym(L(n)*alpha_U(n)))/gamma(sym(L(n)*alpha_U(n))))...
    * v.^(L(n)*alpha_U(n)-1) .* exp( -beta_U(n)*v );                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_V_n = @(v,n) gammainc(beta_U(n)*double(v),L(n)*alpha_U(n),'lower');

% Working on T
% The 1st moment of T
E_T1 = 0;
for n = 1:N_RIS
    for l = 1:L(n)
        E_T1 = E_T1 + E_U_nl_k(1,n);
    end
end

% The 2nd  moment of T
E_T2 = 0;
for n = 1:N_RIS
    tmpA = 0;
    for l = 1:L(n)
        tmpA = tmpA + E_U_nl_k(1,n);
    end
    for ii = n+1:N_RIS
        tmpB = 0;
        for l = 1:L(ii)
            tmpB = tmpB + E_U_nl_k(1,ii);
        end
        E_T2 = E_T2 + 2 * tmpA * tmpB;
    end
end

for n = 1:N_RIS
    tmpC = 0;
    for l = 1:L(n)
        tmpC = tmpC + E_U_nl_k(2,n);
    end
    tmpD = 0;
    for l = 1:L(n)
        for v = (l+1):L(n)
            tmpD = tmpD + 2 * E_U_nl_k(1,n) * E_U_nl_k(1,n);
        end
    end
    E_T2 = E_T2 + tmpC + tmpD;
end

% Working on Z
% The 1st moment of Z
E_Z1 = E_h0_k(1) + E_T1;

% The 2nd moment of Z
E_Z2 = E_h0_k(2) + E_T2 + 2 * E_h0_k(1) * E_T1;

% Parameter of the approximate Gamma distribution of Z
alpha_Z = E_Z1^2/(E_Z2 - E_Z1^2);
beta_Z = E_Z1/(E_Z2 - E_Z1^2);

F_Z_Gamma = @(z) gammainc(z*beta_Z,alpha_Z,'lower'); % Gamma CDF of Z

f_Z_Gamma = @(z) 1/gamma(alpha_Z)*(beta_Z)^alpha_Z...
    * z.^(alpha_Z-1) .* exp( -z*beta_Z ); % Gamma PDF of Z

F_Z2_Gamma = @(z) F_Z_Gamma(sqrt(z)); % CDF of Z^2

%-----Asymptotic of Z^2
F_Z_Gamma_asymp = @(z) (z*beta_Z)^alpha_Z/gamma(alpha_Z+1);
F_Z2_Gamma_asymp= @(z) F_Z_Gamma_asymp(sqrt(z));
%=====

%% code_for_OST

% Working on M_V ( max V_n )
F_M_V = @(x) 1; % CDF of V_M
for k = 1:N_RIS
    F_M_V = @(x) F_M_V(x) .* F_V_n(x,k);
end

M = 100; % number of steps in M-staircase approximation
F_R = @(r) 0; % CDF of R
for m = 1:M
    F_R = @(r) F_R(r)...
        + (F_h0(m/M*r) - F_h0((m-1)/M*r)) .* F_M_V((M-m+1)/M*r);
end
F_R2_Gamma = @(r) F_R(sqrt(r)); % CDF of R^2, e2e SNR of the OST scheme

% %-----Asymptotic of R^2
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
% %=====

%% Test CCT Moment-Matching

%-----CDF of CCT Scheme with Gamma Distribution
figure;
% [y,x] = ecdf(Z2_cct); hold on;
[y,x] = ecdf(Z_cct); hold on;
xs = linspace(0,4,100);

plot(x,y); hold on;
% plot(xs,F_Z2_Gamma(xs),'.','markersize',10); hold on;
plot(xs,F_Z_Gamma(xs),'.','markersize',10); hold on;

xlabel('x')
% ylabel('CDF of Z^2')
ylabel('CDF of Z')
legend('True Distribution','Approximate Distribution',...
    'location','se');

% x0 = 100; y0 = 100; width = 300; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====DONE

%-----PDF of CCT Scheme with Gamma Distribution
figure;
% non-symbolic PDF of Z and Z^2
% f_Z_Gamma = @(z) 1/gamma(alpha_Z)*(beta_Z)^alpha_Z...
%     * z.^(alpha_Z-1) .* exp( -z*beta_Z );
% f_Z2_Gamma = @(z) 1./(2*sqrt(z)) .* f_Z_Gamma(sqrt(z));

number_of_bins = 30;
histogram(Z_cct,number_of_bins,'normalization','pdf'); hold on;

%---symbolic PDF
syms sbl_az sbl_bz sbl_z
symbolic_f_Z_Gamma(sbl_az,sbl_bz,sbl_z) = 1/gamma(sbl_az)*(sbl_bz)^sbl_az ...
    * sbl_z^(sbl_az-1) * exp( - sbl_z*sbl_bz );
plot(xs,double(vpa(symbolic_f_Z_Gamma(alpha_Z,beta_Z,xs))),'linewidth',1); hold on;
%====done

xlabel('x')
ylabel('PDF of Z')
legend('True Distribution', 'Approximate Distribution',...
    'location','ne');

% x0 = 100; y0 = 100; width = 300; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====DONE

%% OST Test Moment-Matching

%-----CDF of R^2
figure;
[y,x] = ecdf(R2_ost); hold on;
xs = linspace(0,10,30);
plot(x,y); hold on;
plot(xs,F_R2_Gamma(xs),'.','markersize',10); hold on;
xlabel('x')
ylabel('CDF of R^2')
legend('True Distribution (from Simulation)','Approximate Distribution',...
    'location','best');
%=====done

%-----PDF of R^2
figure;
f_M_V = @(x) 0;
for n = 1:N_RIS
    func_tmp = @(x) 1;
    for t = 1:N_RIS
        if (n ~= t)
            func_tmp = @(x) func_tmp(x) .* F_V_n(x,t);
        end
    end
    f_M_V = @(x) f_M_V(x) + f_V_n(x,n) .* func_tmp(x);
end
f_R = @(r) 0; M = 100;
for m = 1:M
    f_R = @(r) f_R(r)...
        + ((m/M)*f_h0(m/M*r) - ((m-1)/M)*f_h0((m-1)/M*r))...
        .*F_M_V((M-m+1)/M*r)...
        + (F_h0(m/M*r)- F_h0((m-1)/M*r))...
        .*((M-m+1)/M).*f_M_V((M-m+1)/M*r);
end
f_R2 = @(r) 1./(2*sqrt(r)).*f_R(sqrt(r));

xs = linspace(0,2,100);
number_of_bins = 100;
histogram(R2_ost,number_of_bins,'normalization','pdf'); hold on;
plot(xs,double(vpa(f_R2(sym(xs)))),'linewidth',1.5); hold on;
xlabel('x')
ylabel('PDF of R^2')
legend('True Distribution (from Simulation)','Approximate Distribution',...
    'location', 'best');
%=====done

%% Outage Probability

%-----OP
for isnr = 1:length(snrdB)
    snr = 10^(snrdB(isnr)/10);
    % CCT scheme
    P_out_CCT_Gamma_sim(isnr) = mean(snr*Z2_cct < snr_th);
    P_out_CCT_Gamma_ana(isnr) = F_Z2_Gamma(snr_th/snr);
    P_out_CCT_Gamma_asymp(isnr) = F_Z2_Gamma_asymp(snr_th/snr);
    % OST scheme
    P_out_OST_Gamma_sim(isnr) = mean(snr*R2_ost < snr_th);
    P_out_OST_Gamma_ana(isnr) = F_R2_Gamma(snr_th/snr);
    fprintf('OP, SNR = %d \n', snrdB(isnr));
end

figure;
semilogy(snrdB,P_out_CCT_Gamma_sim,'r^:'); hold on;
semilogy(snrdB,P_out_CCT_Gamma_ana,'r-'); hold on;
semilogy(snrdB,P_out_OST_Gamma_sim,'bo:'); hold on;
semilogy(snrdB,P_out_OST_Gamma_ana,'b-'); hold on;
xlabel('Transmit SNR (dB)');
ylabel('Outage Probability');
legend('CCT (sim.)', 'CCT (ana. with Gamma)',...
    'OST (sim.)', 'OST (ana. with Gamma)',...
    'Location','Best');
axis([-Inf Inf 10^(-5) 10^(0)]);

% x0 = 100; y0 = 100; width = 400; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====

%% Ergodic Capacity
% Note: no EC of OST scheme with Gamma distribution

%-----EC of CCT scheme with Gamma distribution
% input: az = alpha_Z, bz = beta_Z, rho = snr

syms az bz rho
EC_cct_Gamma(az,bz,rho) = 1/gamma(az)/log(2)*2^(az-1)/sqrt(pi)...
    * meijerG( [0],[1/2,1],[az/2,az/2+1/2,0,0,1/2],[], (bz/2)^2/rho);
% input: az = alpha_Z, bz = beta_Z, rho = snr

for isnr = 1:length(snrdB)
    snr = 10^(snrdB(isnr)/10);
    EC_CCT_sim(isnr) = mean(log2(1+snr*Z2_cct));
    %_
    EC_CCT_ana(isnr) = double(vpa(EC_cct_Gamma(alpha_Z,beta_Z,snr)));
    fprintf('Gamma, Ergodic Ca, SNR = %d \n', snrdB(isnr));
end

figure;
plot(snrdB,EC_CCT_sim,'go'); hold on;
plot(snrdB,EC_CCT_ana,'g-'); hold on;
xlabel('Transmit SNR (dB)');
ylabel('Ergodic Capacity (bits/sec/Hz');
legend('CCT (sim.)','CCT (ana. with Gamma)',...
    'Location','Best');

% x0 = 100; y0 = 100; width = 400; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====

%% Network Topology

%-----Topo
fig_topo = figure;
scatter(x_source, y_source, 'bo', 'filled')
hold on
scatter(x_des, y_des, 'go', 'filled')
hold on
scatter(x_ris, y_ris, 'rs', 'filled')
hold on

for kk = 1:N_RIS
    text(x_ris(kk), y_ris(kk), num2str(kk));
    hold on
end

axis([x_area_min x_area_max y_area_min y_area_max])
legend('Source', 'Destination', 'RISs')

x0 = 100; y0 = 100; width = 400; height = 200;
set(gcf,'Position', [x0, y0, width, height]); % plot size
set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====DONE

save('data_Gamma_30_25_20.mat', ...
    'P_out_CCT_Gamma_sim', 'P_out_CCT_Gamma_ana', ...
    'P_out_OST_Gamma_sim', 'P_out_OST_Gamma_ana',...
    'EC_CCT_sim', 'EC_CCT_ana')