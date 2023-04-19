function [vec_fit,error_fit] = add_Gaussians(vec_st,n_FD,Max_iter,num_blocks,x_data,options)
%vec_st is the ongoing fit to the distribution function n_FD
%and we would like to add num_blocks new Gaussians to this list
global amp_max
y_data = n_FD(x_data) - poly_Gauss_approx(vec_st,x_data);

if ~exist(['TMP'], 'dir')
    mkdir(['TMP'])
end

if isempty(vec_st) ~= 1
    [amps_st,stds_st] = convert_from_vec(vec_st);
    AMP_amps = 2*min(abs(amps_st));
    AMP_gamma = 2*min(abs(stds_st));
    if AMP_amps < 1e-5
        AMP_amps = 1;
        AMP_gamma = 1;
    end
else 
    AMP_amps = 1;
    AMP_gamma = 1;
end

amps_vec_ind = 0;
stds_vec_ind = 2*AMP_gamma*sqrt(rand(num_blocks,1) + 1i*rand(num_blocks,1));
vec_ind = convert_to_vec(amps_vec_ind,stds_vec_ind);

error_ind = get_error_fit(merge_vecs(vec_ind,vec_st),x_data,y_data);
writematrix(vec_ind,['TMP/vec_ind_' num2str(0) '.txt'])
writematrix(error_ind,['TMP/error_ind_' num2str(0) '.txt'])

for ind = 1:Max_iter
    amps_vec_ind = (rand(num_blocks,1).*exp(1i*rand(num_blocks,1)*2*pi) - 0.5)*2*AMP_amps;
    stds_vec_ind = 2*AMP_gamma*sqrt(rand(num_blocks,1) + 1i*rand(num_blocks,1));
    vec_ind = convert_to_vec(amps_vec_ind,stds_vec_ind);
    ub = [amp_max*ones(size(amps_vec_ind));amp_max*ones(size(amps_vec_ind));2*ones(size(stds_vec_ind));2*ones(size(stds_vec_ind))];
    lb = [-amp_max*ones(size(amps_vec_ind));-amp_max*ones(size(amps_vec_ind));zeros(size(stds_vec_ind));zeros(size(stds_vec_ind))];
    vec_ind = lsqcurvefit(@(x,t)poly_Gauss_approx(x,t),vec_ind,x_data,y_data,lb,ub,options);
    error_ind = get_error_fit(merge_vecs(vec_ind,vec_st),x_data,y_data);
    disp(['Ongoing error = ' num2str(error_ind)])
    writematrix(vec_ind,['TMP/vec_ind_' num2str(ind) '.txt'])
    writematrix(error_ind,['TMP/error_ind_' num2str(ind) '.txt'])
end

error_best = 100;
for ind = 0:Max_iter
    error_ind = load(['TMP/error_ind_' num2str(ind) '.txt']);
    if error_ind < error_best
        error_best = error_ind;error_fit = error_ind;
        vec_fit = load(['TMP/vec_ind_' num2str(ind) '.txt']);
    end
end
vec_fit = post_process(vec_fit);

rmdir 'TMP' s

end