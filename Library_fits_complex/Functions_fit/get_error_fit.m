function error_fit = get_error_fit(vec,x_data,y_data)

error_fit = max(abs(poly_Gauss_approx(vec,x_data) - y_data));

end