clear;close;clc;

set(0,'DefaultAxesFontSize',24,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontSize',24,'DefaultTextFontName','Arial');
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1e7,...
    'FunctionTolerance',1e-10,'MaxIterations',1e6,'StepTolerance',1e-10,...
    'SpecifyObjectiveGradient',true);

addpath('Functions_fit/')

global amp_max
amp_max = 10;

zF = 10^5;
max_num_funcs = 50;
Error_list = zeros(max_num_funcs,1);

n_FD = @(x) 1 ./ (1/zF*exp(x.^2/2) + 1);
x_data = linspace(0,100,1e4).';
y_data = n_FD(x_data);


for ind = 1:max_num_funcs
    Error_list(ind) = load(['DATA_' num2str(zF)  '/error_best_' num2str(ind) '.txt']);
end


figure(1)
axesH = axes;
%axesH.XAxis.MinorTick       = 'on';
%axesH.YAxis.MinorTick       = 'on';
hold on
plot(1:1:max_num_funcs,log10( Error_list),'-d','LineWidth',1.5,'MarkerSize',14,'Color',[158,1,66]/255,'MarkerFaceColor',[190,186,218]/255)
legend({'z1','z2','z3'})
xlabel('Number of Gaussians $M_0$','Interpreter','latex')
ylabel('Error in $n_p$','Interpreter','latex')
box on
xlim([1 50])
ylim([-8.5 0])

ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.01];
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 2;
% ax.FontWeight = 'normal';
set(gca, 'FontName', 'Arial')

% %% improving the best
% 
% num_funcs = 21;
% vec = load(['DATA_' num2str(zF)  '/vec_best_21.txt']);
% 
% [amps_st,stds_st] = convert_from_vec(vec);
% stds_st = stds_st.*sqrt(log(zF)/0.5);
% vec = convert_to_vec(amps_st,stds_st);
% y_fit = poly_Gauss_approx(vec,x_data);
% 
% disp(['Error in the beginning = ' num2str(max(abs( y_fit - y_data)))])
% 
% figure(2)
% axesH = axes;
% hold on
% 
% plot(x_data,y_data,'-r','LineWidth',2)
% plot(x_data,y_fit,'--b','LineWidth',3)
% 
% box on
% xlim([0 100])
% 
% ax = gca;
% ax.XColor = 'k';
% ax.YColor = 'k';
% ax.TickLength = [0.015 0.01];
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = 2;
% % ax.FontWeight = 'normal';
% set(gca, 'FontName', 'Arial')
% 
% 
% %% massage solution
% vec_imp = massage_sol(vec,n_FD,x_data,options);
% 
% disp(['Error upon massaging = ' num2str(max(abs( poly_Gauss_approx(vec_imp,x_data) - y_data)))])
% 
% %%
% [amps,stds] = convert_from_vec(vec);
% amps_abs = abs(amps);
% 
% vec_matr = zeros((num_funcs - 1)*4,num_funcs);
% error_list = zeros(num_funcs,1);
% 
% tic
% for ind = 1:num_funcs
%     [vec_current,error_current] = remove_Gaussian(vec_imp,ind,n_FD,x_data,options);
%     vec_matr(:,ind) = vec_current;
%     error_list(ind) = error_current;
%     disp(['Current error = ' num2str(error_current)])
% end
% toc
% 
% writematrix(vec_matr,['DATA_post_processing/vec_matr_' num2str(num_funcs)])
% writematrix(error_list,['DATA_post_processing/error_list_' num2str(num_funcs)])
% 
% figure(3)
% axesH = axes;
% hold on
% 
% plot(error_list(:,1),'-x','LineWidth',2)
% 
% 
% box on
% xlim([0 100])
% 
% ax = gca;
% ax.XColor = 'k';
% ax.YColor = 'k';
% ax.TickLength = [0.015 0.01];
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = 2;
% % ax.FontWeight = 'normal';
% set(gca, 'FontName', 'Arial')

% num_funcs = 21;
% error_list = load(['DATA_post_processing/error_list_' num2str(num_funcs) '.txt']);
% vec_matr = load(['DATA_post_processing/vec_matr_' num2str(num_funcs) '.txt']);
% ind_parent = 0;
% target_error = 1e-7;
% for ind = 1:num_funcs
%     if error_list(ind,1) < target_error
%         ind_parent = ind_parent + 1;
%         writematrix(vec_matr(:,ind),['DATA_post_processing/vec_' num2str(num_funcs - 1) '_' num2str(ind_parent)])
%     end
% end

%% Post processing

% 
% num_funcs = 20;
% num_parent = 0;
% for ind = 1:50
%     if exist(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind)  '.txt'], 'file')
%         num_parent = num_parent + 1;
%     end
% end
% disp(['Number of parents = ' num2str(num_parent)]);
% 
% num_offspring = 0;
% 
% for ind_parent = 1:num_parent
%     vec = load(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind_parent)  '.txt']);
%     for ind = 1:num_funcs
%         [vec_current,error_current] = remove_Gaussian(vec,ind,n_FD,x_data,options);
%         if error_current < target_error
%             num_offspring = num_offspring + 1;
%             writematrix(vec_current,['DATA_post_processing/vec_' num2str(num_funcs - 1) '_' num2str(num_offspring) '_' num2str(ind_parent)])
%         end
%     end
% end
% 
% disp(['Offspring number = ' num2str(num_offspring)])
% disp(['Completed the reduction from ' num2str(num_funcs) ' to ' num2str(num_funcs - 1)])
% 

% %%
% num_funcs = 19;
% num_parent = 0;
% vec_best = zeros(4*num_funcs,1);
% error_best = 100;
% for ind = 1:50
%     for sub_ind = 1:50
%         if exist(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind) '_' num2str(sub_ind) '.txt'], 'file')
%             num_parent = num_parent + 1;
%             vec = load(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind) '_' num2str(sub_ind) '.txt']);
%             error_fit = get_error_fit(vec,x_data,n_FD(x_data))
%             if error_fit < error_best
%                 error_best = error_fit;
%                 vec_best = vec;
%             end
%         end
%     end
% end
% disp(['Number of parents = ' num2str(num_parent)]);
% disp(['Ongoing error = ' num2str(error_best)]);
% 
% num_offspring = 0;
% 
% vec = vec_best;
% num_funcs = 19;
% target_error = 1e-7;
% for ind = 1:num_funcs
%     [vec_current,error_current] = remove_Gaussian(vec,ind,n_FD,x_data,options);
%     if error_current < target_error
%         num_offspring = num_offspring + 1;
%         writematrix(vec_current,['DATA_post_processing/vec_' num2str(num_funcs - 1) '_' num2str(num_offspring)])
%     end
% end
% 
% disp(['Offspring number = ' num2str(num_offspring)])
% disp(['Completed the reduction from ' num2str(num_funcs) ' to ' num2str(num_funcs - 1)])
% 

%%
% 
% 
% num_funcs = 18;
% num_parent = 0;
% vec_best = zeros(4*num_funcs,1);
% error_best = 100;
% for ind = 1:50
%     if exist(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind) '.txt'], 'file')
%         num_parent = num_parent + 1;
%         vec = load(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind) '.txt']);
%         error_fit = get_error_fit(vec,x_data,n_FD(x_data));
%         if error_fit < error_best
%             error_best = error_fit;
%             vec_best = vec;
%         end
%     end
% end
% disp(['Number of parents = ' num2str(num_parent)]);
% disp(['Ongoing error = ' num2str(error_best)]);
% 
% num_offspring = 0;
% 
% vec = vec_best;
% num_funcs = 18;
% target_error = 1e-7;
% for ind = 1:num_funcs
%     [vec_current,error_current] = remove_Gaussian(vec,ind,n_FD,x_data,options);
%     if error_current < target_error
%         num_offspring = num_offspring + 1;
%         writematrix(vec_current,['DATA_post_processing/vec_' num2str(num_funcs - 1) '_' num2str(num_offspring)])
%     end
% end
% 
% disp(['Offspring number = ' num2str(num_offspring)])
% disp(['Completed the reduction from ' num2str(num_funcs) ' to ' num2str(num_funcs - 1)])




num_funcs = 17;
num_parent = 0;
vec_best = zeros(4*num_funcs,1);
error_best = 100;
for ind = 1:50
    if exist(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind) '.txt'], 'file')
        num_parent = num_parent + 1;
        vec = load(['DATA_post_processing/vec_' num2str(num_funcs) '_' num2str(ind) '.txt']);
        error_fit = get_error_fit(vec,x_data,n_FD(x_data));
        if error_fit < error_best
            error_best = error_fit;
            vec_best = vec;
        end
    end
end
disp(['Number of parents = ' num2str(num_parent)]);
disp(['Ongoing error = ' num2str(error_best)]);

num_offspring = 0;

vec = vec_best;
num_funcs = 17;
target_error = 1e-7;
for ind = 1:num_funcs
    [vec_current,error_current] = remove_Gaussian(vec,ind,n_FD,x_data,options);
    if error_current < target_error
        num_offspring = num_offspring + 1;
        writematrix(vec_current,['DATA_post_processing/vec_' num2str(num_funcs - 1) '_' num2str(num_offspring)])
    end
end

disp(['Offspring number = ' num2str(num_offspring)])
disp(['Completed the reduction from ' num2str(num_funcs) ' to ' num2str(num_funcs - 1)])




