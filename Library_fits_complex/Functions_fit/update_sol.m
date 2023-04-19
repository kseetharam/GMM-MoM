function vec_fit = update_sol(vec_st,n_FD,Max_iter,x_data,options)
num_blocks = 1;

%Step #1: determine the weakest harmonic
Num_Gaussians_st = length(vec_st)/4;
Error_list = zeros(Num_Gaussians_st,1);
[amps_st,stds_st] = convert_from_vec(vec_st);

for num = 1:Num_Gaussians_st
    amps = amps_st;stds = stds_st;
    amps(num) = [];stds(num) = [];
    Error_list(num) = get_error_fit(convert_to_vec(amps,stds),x_data,n_FD(x_data));
end

[~,index_min] = min(Error_list);
amps_st(index_min) = [];
stds_st(index_min) = [];

%Step #2: massage the ongoing solution plus add two new Gaussians
vec_st = convert_to_vec(amps_st,stds_st);
vec_st = massage_sol(vec_st,n_FD,x_data,options);
vec = add_Gaussians(vec_st,n_FD,Max_iter,num_blocks,x_data,options);

%Step #3: return output
vec_fit = merge_vecs(vec,vec_st);
vec_fit = massage_sol(vec_fit,n_FD,x_data,options);

end