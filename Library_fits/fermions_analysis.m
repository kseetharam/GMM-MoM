clear;close;clc;

addpath('Functions/')

numpts_p = 1001;
xdata = linspace(0,6,numpts_p).';
ydata = 1./(exp(xdata.^2/2) + 1);%this corresponds to zF = 1

num_funcs = 9;
vec = load(['vec_fermi_' num2str(num_funcs) '.txt']);

y_fit = poly_Gauss_approx(vec,xdata);

% figure(1)
% hold on
% plot(xdata,ydata,'-r','LineWidth',2)
% plot(xdata,y_fit,'--b','LineWidth',3)


%% testing on smaller z:

z0 = 0.75;
mF = 0.1;
ydata_z0 = 1./(1/z0*exp(xdata.^2/2/mF) + 1);

amps = vec(1:num_funcs);stds = vec(num_funcs + 1:end);
amps = amps.*((z0).^(1./stds.^2));
vec_z0 = [amps;sqrt(mF)*stds];
y_z0 = poly_Gauss_approx(vec_z0,xdata);

max(abs(y_z0 - ydata_z0))


figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
plot(xdata,ydata_z0,'-r','LineWidth',2)
plot(xdata,y_z0,'--b','LineWidth',3)

