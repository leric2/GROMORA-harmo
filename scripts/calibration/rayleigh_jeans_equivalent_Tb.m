function TbRJE = rayleigh_jeans_equivalent_Tb(calibrationTool, Temperature, frequency)
%==========================================================================
% NAME      | rayleigh_jeans_equivalent_Tb.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2021
%           |
% ABSTRACT  | Function computing the Rayleigh-Jeans equivalent brightness 
%           | temperature. It corresponds to the scaling of the intensity
%           | with unit of temperature. This function convert a Planck or 
%           | physical temperature into a R-J equivalent one (always
%           | lower).
%           |
% ARGUMENTS | INPUTS:   1. calibrationTool: 
%           |               - h, Planck's constant
%           |               - kb, Boltzmann’s constant
%           | 
%           |           2. Temperature: 
%           |               Physical or Planck brightness temperature
%           | 
%           |           3. frequency
%           |
%           | OUTPUTS:  1. TbRJE: Rayleigh-Jeans equivalent brightness 
%           |           temperature. Depending on the inputs, it can be a 
%           |           vector.
%==========================================================================
TbRJE = (calibrationTool.h*frequency ./calibrationTool.kb) ./ (exp(calibrationTool.h * frequency ./ (calibrationTool.kb * Temperature))-1);

end