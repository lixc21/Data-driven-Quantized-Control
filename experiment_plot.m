%% NOTICE FOR READERS
% This code is for paper
% "Data-driven Quantized Control of Partially
% Unknown Linear Systems with Noises"
% 
% Writen by Xingchen Li
% lixc21@mails.tsinghua.edu.cn
% Last modification at 2021-12-07
%
% You should have YALMIP with MOSEK

% Clear variables
clear

%% Plot percentage varies with w_max
% Data preparation
load data_noise_low
w_list=exp((0:.01:1)*-log(0.01))*0.01;
% Draw a semilogarithmic graph
semilogx(w_list(1:5:end),sum(d(1:5:end,:)'>1e-7)/1000,'linewidth',2);
% Set figure properties
xlabel('ω_{max}');
ylabel('percentage')
grid on
set(gca,'FontSize',15);
set(gca,'linewidth',1.5);
xticks([0.01 0.03 0.1 0.3 1])
% Save figure as PDF
saveas(gcf,'noise_exist.pdf')
close all

%% Plot average delta^2 varies with w_max
% Data preparation: integrate two experiments
load data_noise_low
d_low=d;
load data_noise_high
w_list_high=exp((0:.025:1)*(log(1e-2)-log(1e-6)))*1e-6;
w_list_all=[w_list_high(1:end-1) w_list(1:5:65)];
d_all=[d(1:end-1,:); d_low(1:5:65,:)];
% Draw a semilogarithmic graph
semilogx(w_list_all,sum(d_all')./sum(d_all'>1e-7),'linewidth',2);
% Set figure properties
xlabel('ω_{max}');
ylabel('average δ^2')
set(gca,'FontSize',15);
set(gca,'linewidth',1.5);
grid on
xticks([0.000001,0.00001,0.0001,0.001,0.01,0.1])
yticks(0:.05:0.35)
xlim([1e-6 0.15])
% Save figure as PDF
saveas(gcf,'noise_delta.pdf')
close all

%% Plot percentage varies with zeta
% Data preparation
load data_bound
gamma_list=exp((0:.01:1)*log(50));
% Draw a semilogarithmic graph
semilogx(gamma_list(1:5:end),sum(d(1:5:end,:)'>1e-7)/1000,'linewidth',2);
% Set figure properties
xlabel('\zeta');
ylabel('percentage')
set(gca,'FontSize',15);
set(gca,'linewidth',1.5);
grid on
xlim([1 50])
xticks([1 2 3 5 7 10 20 30 50])
% Save figure as PDF
saveas(gcf,'bound_exist.pdf')
close all

%% Plot average delta^2 varies with zeta
% Draw a semilogarithmic graph
semilogx(gamma_list(1:5:91),sum(d(1:5:91,:)')./sum(d(1:5:91,:)'>1e-7),'linewidth',2);
% Set figure properties
xlabel('\zeta');
ylabel('average δ^2')
set(gca,'FontSize',15);
set(gca,'linewidth',1.5);
grid on
xticks([1 2 3 5 7 10 20 32])
xlim([1 32])
% Save figure as PDF
saveas(gcf,'bound_delta.pdf')
close all





