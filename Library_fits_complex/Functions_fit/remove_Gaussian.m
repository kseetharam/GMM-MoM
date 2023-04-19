function [vec_fit,error_fit] = remove_Gaussian(vec,ind,n_FD,x_data,options)

[amps,stds] = convert_from_vec(vec);
amps(ind) = [];
stds(ind) = [];

vec_fit = massage_sol(convert_to_vec(amps,stds),n_FD,x_data,options);
error_fit = get_error_fit(vec_fit,x_data,n_FD(x_data));

end