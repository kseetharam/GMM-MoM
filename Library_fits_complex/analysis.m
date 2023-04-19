clear;close;clc;

set(0,'DefaultAxesFontSize',24,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontSize',24,'DefaultTextFontName','Arial');

addpath(['Functions_fit/'])

numpts_p = 1001;
xdata = linspace(0,6,numpts_p).';

zFmax = 10^5;
vec_f = load(['DATA_' num2str(zFmax)  '/vec_best_21.txt']);
[amps_fermi,stds_fermi] = convert_from_vec(vec_f);
stds_fermi = stds_fermi.*sqrt(log(zFmax));%note this additional factor

mF = 0.1;
ydata_zFmax = 1./(1/zFmax*exp(xdata.^2/mF) + 1);
vec_zFmax = convert_to_vec(amps_fermi,sqrt(mF)*stds_fermi);
y_zFmax = poly_Gauss_approx(vec_zFmax,xdata);


figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
plot(xdata,ydata_zFmax,'-r','LineWidth',2)
plot(xdata,y_zFmax,'--b','LineWidth',3)



zF = 10^1;%choose your own fugacity that is smaller than zF_max
ydata_zF = 1./(1/zF*exp(xdata.^2/mF) + 1);
vec_zF = convert_to_vec(amps_fermi.*((zF/zFmax).^(1./2./stds_fermi.^2)),sqrt(mF)*stds_fermi);
y_zF = poly_Gauss_approx(vec_zF,xdata);

max(abs(y_zF - ydata_zF))


figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
plot(xdata,ydata_zF,'-r','LineWidth',2)
plot(xdata,y_zF,'--b','LineWidth',3)

