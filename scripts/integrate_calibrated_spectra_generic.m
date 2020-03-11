function level1b = integrate_calibrated_spectra_generic(retrievalTool,calibratedSpectra)
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
dt=retrievalTool.integrationTime/60;
timeThresh=0:dt:23;

level1b=struct();
level1b.calibration=calibratedSpectra;

for h = 1:length(timeThresh)
    % Finding the spectra during this time stamp:
    indSpectra=find([calibratedSpectra.TOD]>timeThresh(h) & [calibratedSpectra.TOD]<timeThresh(h)+dt);
    
    % Selecting only spectra with:
    
    % transmission > 0.2 in the integration:
    goodSpectra=indSpectra([calibratedSpectra(indSpectra).meanTroposphericTransmittance] > 0.2);
    
    if ~isempty(goodSpectra)  
        % no critical error:
        goodSpectra=goodSpectra(sum(vertcat(calibratedSpectra(goodSpectra).flags)==[1 1 0 0 0 1 1 1 1],2)==9);
    end
    
    if isempty(goodSpectra)
        intSpectra=NaN*ones(1,retrievalTool.numberOfChannels);
        integratedSpectra(h).numberOfAveragedSpectra=length(goodSpectra);
        goodSpectra=indSpectra;
    else
        % Averaging the good spectra together
        intSpectra=mean(vertcat(calibratedSpectra(goodSpectra).Tb),1);
        integratedSpectra(h).numberOfAveragedSpectra=length(goodSpectra);
    end

    integratedSpectra(h).if=calibratedSpectra(indSpectra(1)).if;
    integratedSpectra(h).freq=calibratedSpectra(indSpectra(1)).freq;
    integratedSpectra(h).timeMin=calibratedSpectra(indSpectra(1)).timeMin;
    integratedSpectra(h).calibrationTime=calibratedSpectra(1).calibrationTime;
    integratedSpectra(h).integrationTime=retrievalTool.integrationTime*60;
    integratedSpectra(h).Tb=intSpectra;
    
    integratedSpectra(h).TOD=mean([calibratedSpectra(goodSpectra).TOD]);
    integratedSpectra(h).meanAirTemperature=mean([calibratedSpectra(goodSpectra).meanAirTemperature]);
    integratedSpectra(h).meanRelativeHumidity=mean([calibratedSpectra(goodSpectra).meanRelHumidity]);
    integratedSpectra(h).TWindow=mean([calibratedSpectra(goodSpectra).TWindow]);
    
end
level1b.integration=integratedSpectra;
end