function [vec,abs_amps] = post_process(vec)
%sorts the elements of vec
Num_Gaussians = length(vec)/4;
amps_vec = vec(1:Num_Gaussians) + 1i*vec(Num_Gaussians + 1:2*Num_Gaussians);
stds_vec = vec(2*Num_Gaussians + 1:3*Num_Gaussians) + 1i*vec(3*Num_Gaussians + 1:4*Num_Gaussians);

abs_amps = abs(amps_vec);
[~,inds] = sort(abs_amps);

abs_amps = flip(abs_amps(inds));
amps_vec = flip(amps_vec(inds));
stds_vec = flip(stds_vec(inds));
vec = [real(amps_vec);imag(amps_vec);...
    real(stds_vec);imag(stds_vec)];
% 
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp('Amplitudes and phases')
% disp([abs_amps angle(amps_vec)])
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp('Gammas and their phases')
% disp([abs(stds_vec) angle(stds_vec)])

end