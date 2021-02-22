function TbPlanck = planck_Tb(calibrationTool, intensity, frequency)
%==========================================================================
% NAME      | planck_Tb.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2021
%           |
% ABSTRACT  | Function computing the Planck brightness temperature. It
%           | corresponds to the temperature of a blackbody radiator that
%           | produces the same intensity as the source being observed
%           | see [Janssen]
%           | This is basically the inverse of the Planck function.
%           | 
%           |
% ARGUMENTS | INPUTS:   1. calibrationTool: 
%           |               - h, Planck's constant
%           |               - kb, Boltzmann’s constant
%           |               - lightSpeed
%           |           2. intensity: 
%           |               The intensity of radiation
%           |           3. frequency
%           |
%           | OUTPUTS:  1. TbPlanck: Planck brightness temperature.
%           | Depending on the inputs, it can be a vector.
%           |
%==========================================================================
TbPlanck = ( calibrationTool.h * frequency ./ (calibrationTool.kb) )./...
    log(1 +  (2*calibrationTool.h*frequency.^3) ./ (intensity* calibrationTool.lightSpeed^2) );
end