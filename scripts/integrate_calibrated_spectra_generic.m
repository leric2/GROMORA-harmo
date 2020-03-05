function integratedSpectra = integrate_calibrated_spectra_generic(calibratedSpectra,integrationTime)
%==========================================================================
% NAME          | 
% TYPE          |
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 09.2014
%               |
% ABSTRACT      |
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:
%               |
%               | OUTPUTS:
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
% Threshold for the integration, calibTime has to be in [min]
dt=integrationTime/60;
timeThresh=0:dt:23;

integratedSpectra=struct();

for h = 1:length(timeThresh)
    % Finding the spectra during this time stamp:
    indSpectra=find([calibratedSpectra.TOD]>timeThresh(h) & [calibratedSpectra.TOD]<timeThresh(h)+dt);
    
    % Selecting only spectra with transmission > 0.2 in the integration:
    goodSpectra=indSpectra([calibratedSpectra(indSpectra).meanTroposphericTransmittance] > 0.2);
    
    % Averaging the good spectra together
    intSpectra=mean(vertcat(calibratedSpectra(goodSpectra).Tb),1);
   
    integratedSpectra(h).intTb=intSpectra;
    integratedSpectra(h).numberOfAveragedSpectra=length(goodSpectra);
    integratedSpectra(h).meanTOD=mean([calibratedSpectra(goodSpectra).TOD]);
end