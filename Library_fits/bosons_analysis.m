clear;close;clc;

set(0,'DefaultAxesFontSize',24,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontSize',24,'DefaultTextFontName','Arial');

addpath('Dilute_bosons/','Functions/')
addpath('Functions/')

z1 = 0.9;
z2 = 0.99;

num_funcs_1 = 11;
num_funcs_2 = 17;
vec_1 = load(['vec_bose_' num2str(z1) '_11.txt']);
vec_2 = load(['vec_bose_' num2str(z2) '_17.txt']);


numpts_p = 1001;
xdata = linspace(0,6,numpts_p).';

ydata_1 = 1./(1/z1*exp(xdata.^2/2) - 1);
ydata_2 = 1./(1/z2*exp(xdata.^2/2) - 1);
y_fit_1 = poly_Gauss_approx(vec_1,xdata);
y_fit_2 = poly_Gauss_approx(vec_2,xdata);

max(abs(ydata_1 - y_fit_1))
max(abs(ydata_2 - y_fit_2))

%% any z < zmax

zBmax = z2;num_funcs = num_funcs_2;
vec_bose = vec_2;

zB = 0.99;%choose your own 0<z<zmax
ydata = 1./(1/zB*exp(xdata.^2/2) - 1);

amps_bose = vec_bose(1:num_funcs);stds_bose = vec_bose(num_funcs + 1:end);
vec_zB = [amps_bose.*((zB/zBmax).^(1./stds_bose.^2));stds_bose];
error_zB = max(abs(ydata - poly_Gauss_approx(vec_zB,xdata)))%error of the fit


figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
plot(xdata,ydata,'-k','LineWidth',3)
plot(xdata,poly_Gauss_approx(vec_zB,xdata),'--','LineWidth',5,'Color',[158,1,66]/255)

legend('full','fit')
box on
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.01];
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 2;
% ax.FontWeight = 'normal';
set(gca, 'FontName', 'Arial')




