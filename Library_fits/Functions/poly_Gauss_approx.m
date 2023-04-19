function [vec_fit,J] = poly_Gauss_approx(vec,xdata)
%here we use real harmonics
%J is the Jacobian

num_funcs = length(vec)/2;
amps_vec = vec(1:num_funcs);
stds_vec = vec((num_funcs + 1):2*num_funcs);

vec_fit = zeros(size(xdata));
J_a = zeros(length(xdata),length(amps_vec));
J_g = zeros(length(xdata),length(amps_vec));

for s = 1:num_funcs
        vec_fit = vec_fit + amps_vec(s).*exp(-0.5*xdata.^2/(stds_vec(s)^2));
        J_a(:,s) = exp(-0.5*xdata.^2/(stds_vec(s)^2));
        J_g(:,s) = amps_vec(s).*xdata.^2.*exp(-0.5*xdata.^2/(stds_vec(s)^2))/(stds_vec(s)^3);
end

J = [J_a J_g];

end