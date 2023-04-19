function [vec_fit,error_fit] = add_Gaussian(vec,xdata,ydata,options,flag)
%Adds a single random Gaussian to existing list
%this random Gaussian is being optimized

num_funcs = length(vec)/2; 
amps = vec(1:num_funcs);
stds = vec((num_funcs + 1):2*num_funcs);
error_in = max(abs(ydata - poly_Gauss_approx(vec,xdata)));

if num_funcs == 0 
    AMP_max = 20;
    STD_max = 2;
else
    AMP_max = min(amps);
    STD_max = min(stds);
end

if flag == 1
    lb = [0;0];ub = [20;2];
else
    lb = [-20;0];ub = [20;2];
end

ydata_upd = ydata - poly_Gauss_approx(vec,xdata);
vec_best = [0;1];
error_best = max(abs(ydata_upd - poly_Gauss_approx(vec_best,xdata)));
Max_Iter = 16;%number of random samples generated
for ind = 1:Max_Iter
    vec_in = [rand(1)*AMP_max;rand(1)*STD_max];
    vec_current = lsqcurvefit(@(x,t)poly_Gauss_approx(x,t),vec_in,xdata,ydata_upd,lb,ub,options);
    error_current = max(abs(ydata_upd - poly_Gauss_approx(vec_current,xdata)));
    if error_current < error_best
        vec_best = vec_current;
        error_best = error_current;
    end
end

amps = [amps;vec_best(1)];
stds = [stds;vec_best(2)];

[amps,inds] = sort(amps);
stds = stds(inds);
vec_fit = [flip(amps);flip(stds)];
error_fit = max(abs(ydata - poly_Gauss_approx(vec_fit,xdata)));

disp(['Error in the beginning = ' num2str(error_in) '; error upon adding new Gaussian = ' num2str(error_fit)])

end