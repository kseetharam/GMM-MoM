function [vec_fit,error_fit] = get_Taylor_vec_fit(z0,num_funcs,xdata,ydata,lb,ub,options,flag)
%Returns the fit if the initial state was chosen to be the Taylor expansion

vec_fit = get_vec_in(num_funcs,z0,1,flag);
vec_fit = lsqcurvefit(@(x,t)poly_Gauss_approx(x,t),vec_fit,xdata,ydata,lb,ub,options);
error_fit = max(abs(ydata - poly_Gauss_approx(vec_fit,xdata)));

end