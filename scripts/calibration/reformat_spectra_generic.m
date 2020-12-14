function rawSpectra = reformat_spectra_generic(rawSpectra,logFile,calibrationTool)
%==========================================================================
% NAME          | reformat_spectra_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 12.2020
%               |
% ABSTRACT      | Transforms the rawSpectra from a vector to a matrix             
%               |
% ARGUMENTS     | INPUTS: 1. rawSpectra as vector
%               |         2. logFile
%               |         3. calibrationTool:
%               |           - numberOfChannels
%               |         
%               | OUTPUTS: - rawSpectra as matrix
%               |          
%==========================================================================
% Number of spectra for this day
%n=size(log.t,1);
%n=length(logFile.t);
% reshaping the initial raw vector
try
    %rawSpectra=(reshape(rawSpectra,[calibrationTool.numberOfChannels,n]))';
    rawSpectra=(reshape(rawSpectra,calibrationTool.numberOfChannels,[]))';
catch ME
    error(ME.identifier,'Problem when reformatting the spectra');
    disp(ME.message)
    
end
end

