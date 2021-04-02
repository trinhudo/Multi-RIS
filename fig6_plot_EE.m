clear;
close all;

SE = 1:.5:20; %Target spectral efficiency [b/s/Hz]

figure;
load('data_EE_setL2.mat')
plot(SE, EE_ORA, 'r-', 'linewidth', 1.5); hold on
load('data_EE_setL3.mat')
plot(SE, EE_ORA, 'b-', 'linewidth', 1.5); hold on
load('data_EE_setL4.mat')
plot(SE, EE_ORA, 'g-', 'linewidth', 1.5); hold on

%--------------------------------------------------------------------------

load('data_EE_setL2.mat')
plot(SE, EE_ERA, 'r-.', 'linewidth', 1.5); hold on
load('data_EE_setL3.mat')
plot(SE, EE_ERA, 'b-.', 'linewidth', 1.5); hold on
load('data_EE_setL4.mat')
plot(SE, EE_ERA, 'g-.', 'linewidth', 1.5); hold on

%--------------------------------------------------------------------------

xlabel('Target SE, $R_{\rm th}$ [b/s/Hz]', ...
    'interpreter', 'latex')
ylabel('Energy efficiency [Mbit/Joule]', ...
    'interpreter', 'latex')
legend('Setting $\bf\mathrm{L}_2$',...
    'Setting $\bf\mathrm{L}_3$',...
    'Setting $\bf\mathrm{L}_4$',...
    'Location','ne',...
    'Interpreter', 'Latex');
axis([1 20 0 200])
set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
set(gca,'fontsize',13);