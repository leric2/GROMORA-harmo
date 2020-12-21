function rawSpectra = flip_spectra_gromos(rawSpectra, calibrationTool)
%==========================================================================
% NAME          | flip_spectra_gromos.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 12.2020
%               |
% ABSTRACT      | Flip the spectra on GROMOS on the frequency axis: 
%               | [0:-500MHz and 0:500] --> -500:500MHz
%               | 
% ARGUMENTS     | INPUTS:   1. rawSpectra
%               |
%               | OUTPUTS:  1. rawSpectra
%               |
%==========================================================================
% flip the negative frequencies (the first 16384 channels) to obtain -500:500MHz

midChannel= calibrationTool.flipAroundChannel;
rawSpectra = [fliplr(rawSpectra(:,1:midChannel)) rawSpectra(:,midChannel+1:end)];
end

