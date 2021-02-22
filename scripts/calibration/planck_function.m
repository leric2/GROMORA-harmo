function Brightness = planck_function(calibrationTool, Temperature, frequency)
%==========================================================================
% NAME      | planck_function.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2021
%           |
% ABSTRACT  | Function computing the spectral intensity emitted by a black
%           | body (Planck's law) from Temperature and frequencies.
%           |
% ARGUMENTS | INPUTS:   1. calibrationTool: 
%           |               - h, Planck's constant
%           |               - kb, Boltzmann’s constant
%           |               - lightSpeed
%           |           2. Temperature: 
%           |               The temperature of the blackbody
%           |           3. frequency
%           |
%           | OUTPUTS:  1. Brightness: The spectral intensity or brightness
%           |               of the radiation emitted.
%           |
%==========================================================================
Brightness = ( (2*calibrationTool.h* frequency.^3) / (calibrationTool.lightSpeed^2) ) ./...
    (exp(calibrationTool.h * frequency ./ (calibrationTool.kb * Temperature))-1);
end