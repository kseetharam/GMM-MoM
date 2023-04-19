function [amps_vec,stds_vec] = convert_from_vec(vec)

if length(vec) > 1
    Num_Gaussians = length(vec)/4;
    amps_vec = vec(1:Num_Gaussians) + 1i*vec(Num_Gaussians + 1:2*Num_Gaussians);
    stds_vec = vec(2*Num_Gaussians + 1:3*Num_Gaussians) + 1i*vec(3*Num_Gaussians + 1:4*Num_Gaussians);
else
    amps_vec = [];
    stds_vec = [];
end

end