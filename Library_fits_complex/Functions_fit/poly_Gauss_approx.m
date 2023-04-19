function [my_fit,J] = poly_Gauss_approx(vec,x)

[amps_vec,stds_vec] = convert_from_vec(vec);
my_fit = zeros(size(x));

J_a_star = zeros(length(x),length(amps_vec));
J_g_star = zeros(length(x),length(amps_vec));

for ind = 1:length(amps_vec)
    my_fit = my_fit + 2*real(amps_vec(ind)*exp(-0.5*x.^2./(stds_vec(ind)^2)));
    J_a_star(:,ind) = conj(exp(-0.5*x.^2./(stds_vec(ind)^2)));
    J_g_star(:,ind) = conj(amps_vec(ind)*exp(-0.5*x.^2./(stds_vec(ind)^2)).*x.^2/(stds_vec(ind)^3));
end

J = [2*real(J_a_star) 2*imag(J_a_star) 2*real(J_g_star) 2*imag(J_g_star)];

end