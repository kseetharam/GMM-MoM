clear;close;clc;
format long;

global amp_max
amp_max = 10;
addpath('Functions_fit/')
set(0,'DefaultAxesFontSize',24,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontSize',24,'DefaultTextFontName','Arial');
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1e7,...
    'FunctionTolerance',1e-10,'MaxIterations',1e6,'StepTolerance',1e-10,...
    'SpecifyObjectiveGradient',true);

z = 1e4;
n_FD = @(x) 1 ./ (1/z*exp(x.^2*log(z)) + 1);
x_data = linspace(0,6,1e4).';
error_tol = 1e-8;
Max_iter = 64;
num_blocks = 1;
Num_Gaussians = 50;
Error_list = zeros(Num_Gaussians,1);

if ~exist(['DATA_' num2str(z)], 'dir')
    mkdir(['DATA_' num2str(z)])
end


vec_st = [];
tic
if ~exist(['DATA_' num2str(z) '/vec_best_1.txt'], 'file')
    tic
    [vec,error_fit] = add_Gaussians(vec_st,n_FD,Max_iter,num_blocks,x_data,options);
    toc
    writematrix(vec,['DATA_' num2str(z) '/vec_best_1.txt'])
    writematrix(error_fit,['DATA_' num2str(z) '/error_best_1.txt'])
    Error_list(1) = error_fit;
else
    vec = load(['DATA_' num2str(z) '/vec_best_1.txt']);
    Error_list(1) = load(['DATA_' num2str(z) '/error_best_1.txt']);
end


vec_st_1 = vec;
vec_st = vec;

for num = 2:Num_Gaussians
    if Error_list(num - 1) < error_tol
        break;
    end
    if ~exist(['DATA_' num2str(z) '/vec_best_' num2str(num) '.txt'], 'file')
        vec = add_Gaussians(vec_st,n_FD,Max_iter,num_blocks,x_data,options);
%         %routine #1
%         vec_st = post_process(merge_vecs(vec,vec_st));N_st = length(vec_st)/4;
%         ub = [1e1*ones(N_st,1);1e1*ones(N_st,1);2*ones(N_st,1);2*ones(N_st,1)];
%         lb = [-1e1*ones(N_st,1);-1e1*ones(N_st,1);zeros(N_st,1);zeros(N_st,1)];
%         vec_st = lsqcurvefit(@(x,t)poly_Gauss_approx(x,t),vec_st,x_data,n_FD(x_data),lb,ub,options);
%         %routine #2
%         vec_st = re_optimize(vec,post_process(vec_st),n_FD,x_data,options);
        %routine #3
        error_prev = get_error_fit(vec_st,x_data,n_FD(x_data));
        vec_st_prev = vec_st;
        vec_st = massage_sol(merge_vecs(vec_st,vec),n_FD,x_data,options);
        error_fit = get_error_fit(vec_st,x_data,n_FD(x_data));
        if (error_fit > 0.7*error_prev)
            vec_st = update_sol(vec_st,n_FD,Max_iter,x_data,options);
            error_fit = get_error_fit(vec_st,x_data,n_FD(x_data));
        end
        Error_list(num) = error_fit;
        writematrix(vec_st,['DATA_' num2str(z) '/vec_best_' num2str(num) '.txt'])
        writematrix(error_fit,['DATA_' num2str(z) '/error_best_' num2str(num) '.txt'])
    else
        vec_st = post_process(load(['DATA_' num2str(z) '/vec_best_' num2str(num) '.txt']));
        Error_list(num) = load(['DATA_' num2str(z) '/error_best_' num2str(num) '.txt']);
    end
end
toc


writematrix(Error_list,['DATA_' num2str(z) '/Error_list.txt'])

figure(1)
axesH = axes;
axesH.XAxis.MinorTick       = 'on';
axesH.YAxis.MinorTick       = 'on';
hold on
plot(log10(Error_list),'-d','LineWidth',1.5,'MarkerSize',14,'Color',[44,127,184]/255,'MarkerFaceColor',[255,255,204]/255)
xlabel('number of Gaussian pairs','Interpreter','latex')
ylabel('$\log L$','Interpreter','latex')
% title(['$z = $' num2str(z) '; without optimization'],'Interpreter','latex')
title(['$z = $' num2str(z) '; with optimization'],'Interpreter','latex')
box on
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.01];
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
% ax.FontWeight = 'normal';
set(gca, 'FontName', 'Arial')

figure(2)
axesH = axes;
axesH.XAxis.MinorTick       = 'on';
axesH.YAxis.MinorTick       = 'on';
plot(x_data,n_FD(x_data) - poly_Gauss_approx(vec_st,x_data),'-r')

