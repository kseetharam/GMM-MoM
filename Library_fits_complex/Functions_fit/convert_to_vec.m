function vec = convert_to_vec(amps_vec,stds_vec)

vec = [real(amps_vec);...
    imag(amps_vec);
    real(stds_vec);
    imag(stds_vec)];

end