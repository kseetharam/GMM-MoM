clear;close;clc;

addpath('Functions/')
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1e7,...
    'FunctionTolerance',1e-12,'MaxIterations',1e6,'StepTolerance',1e-9,...
    'SpecifyObjectiveGradient',true);


z0 = 0.3;

numpts_p = 1001;
xdata = linspace(0,6,numpts_p).';
ydata = 1./(1/z0*exp(xdata.^2/2) - 1);


%% Part 1: collect for the initial state corresponding to the Taylor series

flag = 1;
max_num_funcs = 4;
for num_funcs = 1:max_num_funcs
    ub = [25*ones(num_funcs,1);2*ones(num_funcs,1)];%upper boundary
    lb = [0*ones(num_funcs,1);zeros(num_funcs,1)];%lower boundary
    
    if exist(['Dilute_bosons/vec_bose_TS_' num2str(z0) '_' num2str(num_funcs) '.txt'], 'file')
            vec = load(['Dilute_bosons/vec_bose_TS_' num2str(z0) '_' num2str(num_funcs) '.txt']);
            disp(['z = ' num2str(z0) ';Error_fit = ' num2str( max(abs(ydata - poly_Gauss_approx(vec,xdata))) )]);
    else
        [vec,error_fit] = get_Taylor_vec_fit(z0,num_funcs,xdata,ydata,lb,ub,options,flag);
        disp(['Error_fit = ' num2str( error_fit )]);
        writematrix(vec,['Dilute_bosons/vec_bose_TS_' num2str(z0) '_' num2str(num_funcs) '.txt'])
    end
end
disp('BOSONS ARE COMPLETE')

%% Part 2: finding absolute best within several approaches

Error_list = zeros(num_funcs,1);
for num_funcs = 1:max_num_funcs
    ub = [20*ones(num_funcs,1);2*ones(num_funcs,1)];
    lb = [0*ones(num_funcs,1);zeros(num_funcs,1)];
    if num_funcs == 1
        vec = load(['Dilute_bosons/vec_bose_TS_' num2str(z0) '_' '1.txt']);
        writematrix(vec,['Dilute_bosons/vec_bose_' num2str(z0) '_' num2str(num_funcs) '.txt']);
        Error_list(1) = max(abs(ydata - poly_Gauss_approx(vec,xdata)));
    else
        if ~exist(['Dilute_bosons/vec_bose_' num2str(z0) '_' num2str(num_funcs) '.txt'], 'file')%it loads previous best
            vec_prev = load(['Dilute_bosons/vec_bose_' num2str(z0) '_' num2str(num_funcs - 1) '.txt']);
            amps_prev = vec_prev(1:num_funcs - 1);stds_prev = vec_prev(num_funcs:end);     
%         %we will consider four states to be compared:
%         %1) the one obtained from Taylor series (this one has been precomputed)
%         %2) the one obtained by massaging the solution from 1)
%         %3) we add a random Gaussian to vec_prev, then optimize 
%         %4) we add a not so random Gaussian to vec_prev, then optimize
            vec_1 = load(['Dilute_bosons/vec_bose_TS_' num2str(z0) '_' num2str(num_funcs) '.txt']);
            error_1 = max(abs(ydata - poly_Gauss_approx(vec_1,xdata)));
            vec = vec_1;error_best = error_1;
            
            vec_2 = massage_sol(vec_1,xdata,ydata,options,flag);
            error_2 = max(abs(ydata - poly_Gauss_approx(vec_2,xdata)));
            if error_best > error_2
                vec = vec_2;error_best = error_2;
                disp('Massaging TS solution helped')
            end

            vec_3 = add_Gaussian(vec_prev,xdata,ydata,options,flag);
            vec_3 = massage_sol(vec_3,xdata,ydata,options,flag);
            error_3 = max(abs(ydata - poly_Gauss_approx(vec_3,xdata)));
            if error_best > error_3
                vec = vec_3;error_best = error_3;
                disp('New random Gaussian helped')
            end

            vec_4 = [amps_prev;1;stds_prev;1/sqrt(num_funcs)];
            error_4 = max(abs(ydata - poly_Gauss_approx(vec_4,xdata)));
            if error_best > error_4
                vec = vec_4;error_best = error_4;
                disp('New not so random Gaussian helped')
            end
            Error_list(num_funcs) = error_best;
            writematrix(vec,['Dilute_bosons/vec_bose_' num2str(z0) '_' num2str(num_funcs) '.txt'])
            writematrix(Error_list,['Dilute_bosons/Error_list_' num2str(z0) '.txt'])
            disp(['Ongoing error = ' num2str(error_best) '; number of Gaussians = ' num2str(num_funcs)])
        else
            vec = load(['Dilute_bosons/vec_bose_' num2str(z0) '_' num2str(num_funcs) '.txt']);
            Error_list(num_funcs) = max(abs(ydata - poly_Gauss_approx(vec,xdata)));
        end  
    end
end

