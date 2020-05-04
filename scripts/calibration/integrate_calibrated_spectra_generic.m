function level1b = integrate_calibrated_spectra_generic(calibrationTool,level1b)
%==========================================================================
% NAME          | integrate_calibrated_spectra_generic.m
% TYPE          | Function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 03.2020
%               |
% ABSTRACT      | Function for integrating calibrated spectra from an
%               | instrument
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
dt=calibrationTool.integrationTime/60;
timeThresh=0:dt:23;

for h = 1:length(timeThresh)
    % Finding the spectra during this time stamp:
    indSpectra=find([level1b.calibratedSpectra.TOD]>timeThresh(h) & [level1b.calibratedSpectra.TOD]<timeThresh(h)+dt);
    
    % Selecting only spectra with:
    
    % transmission > 0.2 in the integration:
    goodSpectra=indSpectra([level1b.calibratedSpectra(indSpectra).meanTroposphericTransmittance] > 0.2);
    
    if ~isempty(goodSpectra)  
        % no critical error:
        goodSpectra=goodSpectra(sum(vertcat(level1b.calibratedSpectra(goodSpectra).flags)==[1 1 1 1 1 1],2)==6);
    end
    
    if isempty(goodSpectra)
        intSpectra=-9999*ones(1,calibrationTool.numberOfChannels);
        integratedSpectra(h).numberOfAveragedSpectra=-9999;
        goodSpectra=indSpectra;
    else
        % Averaging the good spectra together
        intSpectra=mean(vertcat(level1b.calibratedSpectra(goodSpectra).Tb),1);
        integratedSpectra(h).numberOfAveragedSpectra=length(goodSpectra);
        
        % Summing the number of spectra for hot, cold and antenna:
        integratedSpectra(h).numHotSpectra = sum(level1b.calibratedSpectra(goodSpectra).numHotSpectra);
        integratedSpectra(h).numColdSpectra = sum(level1b.calibratedSpectra(goodSpectra).numColdSpectra);
        integratedSpectra(h).numAntSpectra = sum(level1b.calibratedSpectra(goodSpectra).numAntSpectra);
    end

    integratedSpectra(h).if=level1b.calibratedSpectra(indSpectra(1)).if;
    integratedSpectra(h).freq=level1b.calibratedSpectra(indSpectra(1)).freq;
    integratedSpectra(h).timeMin=level1b.calibratedSpectra(indSpectra(1)).timeMin;
    
    integratedSpectra(h).firstSkyTime=level1b.calibratedSpectra(indSpectra(1)).firstSkyTime;
    integratedSpectra(h).lastSkyTime=level1b.calibratedSpectra(indSpectra(end)).lastSkyTime;
    
    integratedSpectra(h).year=level1b.calibratedSpectra(indSpectra(1)).year;
    integratedSpectra(h).month=level1b.calibratedSpectra(indSpectra(1)).month;
    integratedSpectra(h).day=level1b.calibratedSpectra(indSpectra(1)).day;
    
    integratedSpectra(h).calibrationTime=level1b.calibratedSpectra(1).calibrationTime;
    integratedSpectra(h).integrationTime=calibrationTool.integrationTime*60;
    integratedSpectra(h).Tb=intSpectra;
    integratedSpectra(h).meanAngleAntenna=nanmean([level1b.calibratedSpectra(goodSpectra).meanAngleAntenna]);
    
    integratedSpectra(h).TOD=nanmean([level1b.calibratedSpectra(goodSpectra).TOD]);
    
    
    integratedSpectra(h).meanDatetime=nanmean([level1b.calibratedSpectra(goodSpectra).meanDatetime]);
    
    integratedSpectra(h).meanAirTemperature=nanmean([level1b.calibratedSpectra(goodSpectra).meanAirTemperature]);
    
    if isnan(integratedSpectra(h).meanAirTemperature)
        integratedSpectra(h).meanAirTemperature=-9999;
    end
    
    integratedSpectra(h).meanRelativeHumidity=nanmean([level1b.calibratedSpectra(goodSpectra).meanRelHumidity]);
    if isnan(integratedSpectra(h).meanRelativeHumidity)
        integratedSpectra(h).meanRelativeHumidity=-9999;
    end
    
    integratedSpectra(h).rainAccumulation=nansum([level1b.calibratedSpectra(indSpectra).rainAccumulation]);
    if isnan(integratedSpectra(h).rainAccumulation)
        integratedSpectra(h).rainAccumulation=-9999;
    end
    
    integratedSpectra(h).TSys=nanmean([level1b.calibratedSpectra(goodSpectra).TSys]);
    if isnan(integratedSpectra(h).TSys)
        integratedSpectra(h).TSys=-9999;
    end
    
    integratedSpectra(h).THot=nanmean([level1b.calibratedSpectra(goodSpectra).THot]);
    if isnan(integratedSpectra(h).THot)
        integratedSpectra(h).THot=-9999;
    end
    
    %stdTHot
    %StdTSYS...
    
    integratedSpectra(h).TWindow=nanmean([level1b.calibratedSpectra(goodSpectra).TWindow]);
    integratedSpectra(h).TOut=nanmean([level1b.calibratedSpectra(goodSpectra).TWindow]);
end
level1b.integration=integratedSpectra;
end