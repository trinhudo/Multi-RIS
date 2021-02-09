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

% line 
% x_ris = x_area_min + (x_area_max-x_area_min)*linspace(.1, .9, N_RIS)'; % [num_relay x 1] vector
% y_ris = .7* ones(N_RIS,1)* y_area_max;

% fixed location
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
% MONTE-CARLO SIMULATION
h_0 = random('Naka',m0,Omg0,[1,sim_trials]);

V_n = zeros(N_RIS,sim_trials);
for n = 1:N_RIS
    for l = 1:L(n)
        %
        h_nl = random('Naka',m_h(n),Omgh(n),[1,sim_trials]);
        g_nl = random('Naka',m_g(n),Omgg(n),[1,sim_trials]);
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

%% code_for_CCT. Fit Z2_CCT to Log-Normal

lambda = sqrt( m_g./Omgg.*m_h./Omgh )/kappa_nl;

% Working on h0
% The k-th moment of h0
E_h0_k = @(k) gamma(m0+k/2)/gamma(m0)*(m0/Omg0)^(-k/2);

F_h0 = @(x) gammainc(m0*x.^2/Omg0,m0,'lower');

% Working on U_nl
% The k-moment of U_nl
E_U_nl_k = @(k,n) lambda(n)^(-k)*gamma(m_h(n)+0.5*k)...
    * gamma(m_g(n)+0.5*k) / gamma(m_h(n)) / gamma(m_g(n));

% Parameter of the approximate Gamma distribution of U_nl
alpha_U = @(n) E_U_nl_k(1,n)^2/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);
beta_U = @(n) E_U_nl_k(1,n)/(E_U_nl_k(2,n)-E_U_nl_k(1,n)^2);

% Working on V_n
% The k-moment of V_n
E_V_n_k = @(k,n) gamma(L(n) * alpha_U(n)+k) ...
    / gamma(L(n) * alpha_U(n)) * beta_U(n)^(-k);

% Working on T
% - the 1st moment of T -
E_T1 = 0;
for n = 1:N_RIS
    for l = 1:L(n)
        E_T1 = E_T1 + E_U_nl_k(1,n);
    end
end
% - the 2nd  moment of T -
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

E_Z = E_h0_k(1) + E_T1;
E_Z1 = E_h0_k(1) + E_T1;
E_Z2 = E_h0_k(2) + E_T2 + 2 * E_h0_k(1) * E_T1;

% Fit Z2_CCT to Log-Normal
E_Z4 = 0; % the 2nd moment of Z^2, i.e., E[(Z^2)^2]
for k = 0:4
    %
    E_Tk = 0; % E[T^k] % the k-th momemt of T >> CAN BE USE IN CCT ???
    [cases_T,indT_mat] = nsumk(N_RIS,k);
    for icaseT = 1:cases_T
        indT_arr = indT_mat(icaseT,:);
        %
        tmpT = 1;
        for t = 1:N_RIS
            tmpT = tmpT * E_V_n_k(indT_arr(t),t);
        end
        %
        E_Tk = E_Tk + factorial(k)/prod( factorial(indT_arr) ) * tmpT;
        %
    end
    %
    E_Z4 = E_Z4 + nchoosek(4,k)*E_h0_k(4-k)*E_Tk;
end

nu_cct_LN = log( E_Z2^2/sqrt(E_Z4) ); % nu_Z
zeta2_cct_LN = log( E_Z4/E_Z2^2 ); % zeta_Z^2

nu_cct_new = log(E_Z^2 / sqrt(E_Z2) );
zeta2_cct_new = log( E_Z2 / (E_Z^2) );

F_Z_cct_new = ...
    @(x) 1/2 + 1/2*erf( (log(x)-nu_cct_new)/sqrt(2*zeta2_cct_new) ); % CDF of Z 

F_Z2_cct_LN = ...
    @(x) 1/2 + 1/2*erf( (log(x)-nu_cct_LN)/sqrt(2*zeta2_cct_LN) ); % CDF of Z^2


%% code_for_OST. Fit R and R^2 to Log-Normal

alpha_U_arr = zeros(1,N_RIS);
beta_U_arr = zeros(1,N_RIS);

for n = 1:N_RIS
    alpha_U_arr(n) = alpha_U(n);
    beta_U_arr(n) = beta_U(n);
end
%
chi_t= @(t) beta_U_arr(t)./sum(beta_U_arr);

% Approxiate result, using self-built F_A() function

alpha_U_arr= zeros(1,N_RIS);
beta_U_arr = zeros(1,N_RIS);
for n = 1:N_RIS
    alpha_U_arr(n) = alpha_U(n);
    beta_U_arr(n) = beta_U(n);
end

f_V_n = @(v,n) beta_U(n)^(L(n)*alpha_U(n))/gamma(L(n)*alpha_U(n))...
    * v.^(L(n)*alpha_U(n)-1) .* exp( -beta_U(n)*v );
F_V_n = @(v,n) gammainc(beta_U(n)*v,L(n)*alpha_U(n),'lower');

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

% mu_M_V = zeros(1,4); % the k-th moment of M_V (k = 1,2,3,4)
% for k = 1:4
%     mu_M_V(k) = integral(@(x) x.^k .* f_M_V(x),0,250);
% end

X = sym('X',[1,N_RIS]);
mu_M_V = sym(zeros(1,4)); % the k-th moment of M_V (k = 1,2,3,4)
for k = 1:4
    for n = 1:N_RIS
        %
        Sn = setdiff(1:N_RIS,n);
        %
        tmp = Lauricella_FA(sum(L.*alpha_U_arr)+k,ones(1,N_RIS-1),L(Sn).*alpha_U_arr(Sn)+1,chi_t(Sn));
        %
        if (tmp > 0)
            mu_M_V(k) = mu_M_V(k) + gamma( sum(X)+k ) / gamma(X(n))...
                / prod( gamma(X(Sn)+1) ) * sym(tmp);
        end
    end
    mu_M_V(k) = vpa(subs(sum(beta_U_arr)^(-k) * prod( chi_t(1:N_RIS).^(X) ) * mu_M_V(k),X,L.*alpha_U_arr));
end
mu_M_V = double(mu_M_V);


E_R = 0; % E[R]
for k = 0:1
    if k >= 1
        E_R = E_R + nchoosek(1,k) * E_h0_k(1-k) * mu_M_V(k);
    else
        E_R = E_R + E_h0_k(1);
    end
end

E_R_2 = 0; % E[R^2] by using R^2 expressions, not R
for k = 0:2
    if k >= 1
        E_R_2 = E_R_2 + nchoosek(2,k) * E_h0_k(2-k) * mu_M_V(k);
    else
        E_R_2 = E_R_2 + E_h0_k(2);
    end
end

E_R_4 = 0; % E[(R^2)^2] by using R^2 expressions, not R
for k = 0:4
    if k >= 1
        E_R_4 = E_R_4 + nchoosek(4,k) * E_h0_k(4-k) * mu_M_V(k);
    else
        E_R_4 = E_R_4 + E_h0_k(4);
    end
end

nu_R_ost_LN = log( E_R^2/sqrt(E_R_2) );
zeta2_R_ost_LN = log( E_R_2/E_R^2 ); % zeta^2 of R

nu_R2_ost_LN = log( E_R_2^2/sqrt(E_R_4) );
zeta2_R2_ost_LN = log( E_R_4/E_R_2^2 ); % zeta^2 of R^2

F_R_ost_LN = @(x) 1/2 + ...
    1/2*erf( (log(x)-nu_R_ost_LN)/sqrt(2*zeta2_R_ost_LN) ); % CDF of R

F_R2_ost_LN = @(x) 1/2 + ...
    1/2*erf( (log(x)-nu_R2_ost_LN)/sqrt(2*zeta2_R2_ost_LN) ); % CDF of R^2

%% Test CDF and PDF

% ----------------------- CDF of CCT Scheme with Gamma Distribution -----------------------
figure;
% [y,x] = ecdf(Z2_cct); hold on;
[y,x] = ecdf(Z_cct); hold on;
xs = linspace(0,5,100);

plot(x,y); hold on;
% plot(xs,F_Z2_cct_LN(xs),'.','markersize',10); hold on;
plot(xs,F_Z_cct_new(xs),'.','markersize',10); hold on;


xlabel('x')
% ylabel('CDF of Z^2')
ylabel('CDF of Z')
legend('True CDF', 'LogNormal Approx.',...
    'location','se');

% x0 = 100; y0 = 100; width = 300; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====

% ----------------------- PDF of CCT Scheme with Gamma Distribution -----------------------
figure;

f_Z_cct_new = @(x) 1./(x .* sqrt(2*pi*zeta2_cct_new) )...
    .* exp( -(log(x) - nu_cct_new).^2./(2*zeta2_cct_new) ); % PDF of Z

f_Z2_cct_LN = @(x) 1./(x .* sqrt(2*pi*zeta2_cct_LN) )...
    .* exp( -(log(x) - nu_cct_LN).^2./(2*zeta2_cct_LN) ); % PDF of Z^2

number_of_bins = 30;
% histogram(Z2_cct,number_of_bins,'normalization','pdf'); hold on;
histogram(Z_cct,number_of_bins,'normalization','pdf'); hold on;

plot(xs,double(vpa(f_Z_cct_new(sym(xs)))),'linewidth',1); hold on;
%====done

xlabel('x')
% ylabel('PDF of Z^2')
ylabel('PDF of Z')
legend('True PDF', 'LogNormal Approx.',...
    'location','ne');

% x0 = 100; y0 = 100; width = 300; height = 250;
% set(gcf,'Position', [x0, y0, width, height]); % plot size
% set(gca, 'LooseInset', get(gca, 'TightInset')) % remove plot padding
%=====done

% Test OST Moment-Matching

% --------------------------------- CDF of R -----------------------
figure;
[y,x] = ecdf(R_ost); hold on;
xs = linspace(0,3,30);
plot(x,y); hold on;
plot(xs,F_R_ost_LN(xs),'.','markersize',10); hold on;
xlabel('x')
ylabel('CDF of R')
legend('True CDF','Aprrox. CDF',...
    'location','best');
%=====done

% --------------------------------- PDF of R -----------------------
figure;

f_R_ost_LN = @(x) 1./(x .* sqrt(2*pi*zeta2_R_ost_LN) )...
    .* exp( -(log(x) - nu_R_ost_LN).^2./(2*zeta2_R_ost_LN) ); % CDF of R

xs = linspace(0,2,100);
number_of_bins = 100;
histogram(R_ost,number_of_bins,'normalization','pdf'); hold on;
plot(xs,double(vpa(f_R_ost_LN(sym(xs)))),'linewidth',1.5); hold on;
xlabel('x')
ylabel('PDF of R')
legend('True PDF','Approx. PDF',...
    'location', 'best');
%=====done

% % ---------------------------------------------- CDF of R^2 -----------------------
% figure;
% [y,x] = ecdf(R2_ost); hold on;
% xs = linspace(0,3,30);
% plot(x,y); hold on;
% plot(xs,F_R2_ost_LN(xs),'.','markersize',10); hold on;
% xlabel('x')
% ylabel('CDF of R^2')
% legend('True Distribution (from Simulation)','Approximate Distribution',...
%     'location','best');
% %=====done
% 
% % ---------------------------------------------- PDF of R^2 -----------------------
% figure;
% f_R2_ost_LN = @(x) 1./(x .* sqrt(2*pi*zeta2_ost_LN) )...
%     .* exp( -(log(x) - nu_ost_LN).^2./(2*zeta2_ost_LN) );
% 
% xs = linspace(0,2,100);
% number_of_bins = 100;
% histogram(R2_ost,number_of_bins,'normalization','pdf'); hold on;
% plot(xs,double(vpa(f_R2_ost_LN(sym(xs)))),'linewidth',1.5); hold on;
% xlabel('x')
% ylabel('PDF of R^2')
% legend('True Distribution (from Simulation)','Approximate Distribution',...
%     'location', 'best');
% %=====done


%% OUTAGE PROBABILITY
%-----
for isnr = 1:length(snrdB)
    snr = 10^(snrdB(isnr)/10);
    %
    P_out_CCT_LN_sim(isnr) = mean(snr*Z2_cct < snr_th);
    P_out_CCT_LN_ana(isnr) = F_Z_cct_new(sqrt(snr_th/snr)); % F_Z (sqrt(x))
    %
    P_out_OST_LN_sim(isnr) = mean(snr*R2_ost < snr_th);
    P_out_OST_LN_ana(isnr) = F_R_ost_LN(sqrt(snr_th/snr));
    fprintf('LogNormal, Outage Probability, SNR = %d \n', snrdB(isnr));
end

figure;
semilogy(snrdB,P_out_CCT_LN_sim,'r^:'); hold on;
semilogy(snrdB,P_out_CCT_LN_ana,'r-'); hold on;
semilogy(snrdB,P_out_OST_LN_sim,'bo:'); hold on;
semilogy(snrdB,P_out_OST_LN_ana,'b-'); hold on;
xlabel('Transmit SNR (dB)');
ylabel('Outage Probability');
legend('CCT (Sim)','CCT (Ana-Log-Normal)',...
    'OST (Sim)','OST (Ana-Log-Normal)',...
    'Location','Best');
axis([-Inf Inf 10^(-5.5) 10^(0)]);

% %=====

%% ERGODIC CAPACITY
%-----
a1 = 0.9999964239;
a2 =-0.4998741238;
a3 = 0.3317990258;
a4 =-0.2407338084;
a5 = 0.1676540711;
a6 =-0.0953293897;
a7 = 0.0360884937;
a8 =-0.0064535442;
a_arr= [a1,a2,a3,a4,a5,a6,a7,a8];

syms aXi bXi Xi(aXi,bXi)
Xi(aXi,bXi) = 0;
for k = 1:8
    Xi(aXi,bXi) = Xi(aXi,bXi) + exp(-bXi^2)/2 * a_arr(k)...
        * exp((k/(2*aXi) + bXi)^2) * erfc(k/(2*aXi) + bXi);
end

syms zeta2 nu
EC_LN(zeta2,nu) = ...
    Xi(1/sqrt(2*zeta2),  nu/sqrt(2*zeta2))...
    + Xi(1/sqrt(2*zeta2), -nu/sqrt(2*zeta2))...
    + sqrt(zeta2/2/pi) * exp(-nu^2/(2*zeta2))...
    + nu/2*erfc(-nu/sqrt(2*zeta2));

for isnr = 1:length(snrdB)
    snr = 10^(snrdB(isnr)/10);
    %
    EC_CCT_LN_sim(isnr) = mean(log2(1+snr*Z2_cct));
    EC_CCT_LN_ana(isnr) = double(vpa(EC_LN(zeta2_cct_LN,log(snr) + nu_cct_LN)))/log(2);
    
    EC_OST_LN_sim(isnr) = mean(log2(1+snr*R2_ost));
    EC_OST_LN_ana(isnr) = double(vpa(EC_LN(zeta2_R2_ost_LN,log(snr) + nu_R2_ost_LN)))/log(2);
    fprintf('LogNormal, Ergodic Capacity , SNR = %d \n', snrdB(isnr));
end

figure;
plot(snrdB,EC_CCT_LN_sim,'r*:'); hold on;
plot(snrdB,EC_CCT_LN_ana,'r-'); hold on;
plot(snrdB,EC_OST_LN_sim,'bo:'); hold on;
plot(snrdB,EC_OST_LN_ana,'b-'); hold on;
xlabel('Transmit SNR (dB)');
ylabel('Ergodic Capacity');
legend('CCT (Sim)','CCT (Ana-Log-Normal)',...
    'OST (Sim)','OST (Ana-Log-Normal)',...
    'Location','Best');
%=====

%% Plot topo

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

% save('data_LogNormal_20_25_30.mat', ...
%     'P_out_CCT_LN_sim', 'P_out_CCT_LN_ana', ...
%     'P_out_OST_LN_sim', 'P_out_OST_LN_ana', ...
%     'EC_CCT_LN_sim', 'EC_CCT_LN_ana', ...
%     'EC_OST_LN_sim', 'EC_OST_LN_ana')