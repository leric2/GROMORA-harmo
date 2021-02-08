function B = planck_function(calibrationTool, T, f)

B = ( (2*calibrationTool.h* f.^3) / (calibrationTool.lightSpeed^2) ) ./ (exp(calibrationTool.h * f ./ (calibrationTool.kb * T))-1);

end