function [vec_fit,error_fit] = massage_sol(vec_st,xdata,ydata,options,flag)
%Updates solution

num_funcs = length(vec_st)/2;
error_in = max(abs(ydata - poly_Gauss_approx(vec_st,xdata)));
vec = [];

for ind = 1:num_funcs
    num_funcs_vec = length(vec)/2;
    num_funcs_vec_st = length(vec_st)/2;
    if num_funcs_vec > 0
        amps = vec(1:num_funcs_vec);stds = vec(num_funcs_vec + 1:end);
    else
        amps = [];stds = [];
    end
    amps_st = vec_st(1:num_funcs_vec_st);stds_st = vec_st(num_funcs_vec_st + 1:end);
    
    index = randi([1 num_funcs_vec_st]);
    vec = [amps;amps_st(index);stds;stds_st(index)];
    amps_st(index) = [];stds_st(index) = [];
    vec_st = [amps_st;stds_st];
    num_funcs_vec = length(vec)/2;

    ydata_ind = ydata - poly_Gauss_approx(vec_st,xdata);
    if flag == 1
        ub = [20*ones(num_funcs_vec,1);2*ones(num_funcs_vec,1)];
        lb = [0*ones(num_funcs_vec,1);zeros(num_funcs_vec,1)];
    else
        ub = [20*ones(num_funcs_vec,1);2*ones(num_funcs_vec,1)];
        lb = [-20*ones(num_funcs_vec,1);zeros(num_funcs_vec,1)];
    end

    vec = lsqcurvefit(@(x,t)poly_Gauss_approx(x,t),vec,xdata,ydata_ind,lb,ub,options);
end

num_funcs_vec = length(vec)/2;
amps = vec(1:num_funcs_vec);stds = vec(num_funcs_vec + 1:end);

[amps,inds] = sort(amps);
stds = stds(inds);
vec_fit = [flip(amps);flip(stds)];
error_fit = max(abs(ydata - poly_Gauss_approx(vec_fit,xdata)));

disp(['Error in the beginning = ' num2str(error_in) '; error upon massaging solution = ' num2str(error_fit)])

end