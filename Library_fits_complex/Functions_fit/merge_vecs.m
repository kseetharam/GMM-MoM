function [vec] = merge_vecs(vec_1,vec_2)

[amps_1,stds_1] = convert_from_vec(vec_1);
[amps_2,stds_2] = convert_from_vec(vec_2);

amps = [amps_1;amps_2];stds = [stds_1;stds_2];
vec = convert_to_vec(amps,stds);

end