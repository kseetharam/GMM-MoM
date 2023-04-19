function vec_fit = massage_sol(vec_st,n_FD,x_data,options)
%This function merges vec and vec_st into a single vec_fit, but it does it
%in such a way so that in between it performs a little bit of optimization
global amp_max
Num_Gaussians_st = length(vec_st)/4;
vec = [];

for num = 1:Num_Gaussians_st
    [amps,stds] = convert_from_vec(vec);
    [amps_st,stds_st] = convert_from_vec(vec_st);
    %routine #1
    %amps_st = flip(amps_st); stds_st = flip(stds_st);%this is to avoid oscillations in shallow minima
    %vec = convert_to_vec([amps;amps_st(end)],[stds;stds_st(end)]);

    %routine #2
    Num_st = length(vec_st)/4;
    index = randi([1 Num_st]);
    vec = convert_to_vec([amps;amps_st(index)],[stds;stds_st(index)]);
    Num_Gaussians_vec = length(vec)/4;
    amps_st(index) = [];stds_st(index) = [];
    vec_st = convert_to_vec(amps_st,stds_st);
    y_data = n_FD(x_data) - poly_Gauss_approx(vec_st,x_data);
    ub = [amp_max*ones(Num_Gaussians_vec,1);amp_max*ones(Num_Gaussians_vec,1);2*ones(Num_Gaussians_vec,1);2*ones(Num_Gaussians_vec,1)];
    lb = [-amp_max*ones(Num_Gaussians_vec,1);-amp_max*ones(Num_Gaussians_vec,1);zeros(Num_Gaussians_vec,1);zeros(Num_Gaussians_vec,1)];
    vec = lsqcurvefit(@(x,t)poly_Gauss_approx(x,t),vec,x_data,y_data,lb,ub,options);
end

vec_fit = post_process(vec);
end