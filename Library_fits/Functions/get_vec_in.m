function vec_in = get_vec_in(num_funcs,z,mF,flag)
%Returns initial vector that follows the Taylor series expansion

amplitudes = zeros(num_funcs,1);
gammas = zeros(num_funcs,1);
for ind = 1:num_funcs
    if flag == 1%bosons
        amplitudes(ind) = z^ind;
        gammas(ind) = 1/sqrt(ind);
    else%fermions
        amplitudes(ind) = (-1)^(ind + 1)*z^ind;
        gammas(ind) = sqrt(mF)/sqrt(ind);
    end
end

vec_in = [amplitudes;gammas];

end