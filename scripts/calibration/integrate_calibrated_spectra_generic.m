function integratedSpectra = integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra)
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
dt=hours(calibrationTool.integrationTime/60);
%timeThresh=0:dt:23;
timeThresh = datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day):dt:datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day+1);

for h = 1:length(timeThresh)-1
    % Finding the spectra during this time stamp:
    indSpectra=find([calibratedSpectra.dateTime]>=timeThresh(h) & [calibratedSpectra.dateTime]<timeThresh(h)+dt);
    
    % Selecting only spectra with:
    
    % transmission > 0.2 in the integration:
    if calibrationTool.filterByTransmittance
        goodSpectra=indSpectra([calibratedSpectra(indSpectra).troposphericTransmittance] > 0.2);
    else
        goodSpectra = indSpectra;
    end
    if calibrationTool.filterByFlags
        if ~isempty(goodSpectra)  
            noErrorVect = ones(1,calibrationTool.flagVectorLength);
            % no critical error:
            goodSpectra=goodSpectra(sum(vertcat(calibratedSpectra(goodSpectra).flags)==noErrorVect,2)==calibrationTool.flagVectorLength);
        end
    end
    if isempty(goodSpectra)
        integratedTb=-9999*ones(1,calibrationTool.numberOfChannels);
        integratedStdTb = -9999*ones(1,calibrationTool.numberOfChannels);
        integratedMeanStdTbFromCalMean = -9999;
        integratedSpectra(h).numberOfAveragedSpectra = 0;
        % goodSpectra=indSpectra;
        
        integratedSpectra(h).numHotSpectra = 0;
        integratedSpectra(h).numColdSpectra = 0;
        integratedSpectra(h).numSkySpectra = 0;
        
        meanAngleAT = -9999;
        tod = nanmean([calibratedSpectra(indSpectra).TOD]);
        dateTime=nanmean([calibratedSpectra(indSpectra).dateTime]);
    else
        % Averaging the good spectra together
        integratedTb=nanmean(vertcat(calibratedSpectra(goodSpectra).Tb),1);
        integratedStdTb = nanmean(vertcat(calibratedSpectra(goodSpectra).stdTb),1) / sqrt(length(goodSpectra));
        integratedSpectra(h).numberOfAveragedSpectra=length(goodSpectra);
        
        integratedMeanStdTbFromCalMean = nanmean([calibratedSpectra(goodSpectra).meanStdTb])/sqrt(length(goodSpectra));
        
        % Summing the number of spectra for hot, cold and antenna:
        integratedSpectra(h).numHotSpectra = sum(vertcat(calibratedSpectra(goodSpectra).numHotSpectra));
        integratedSpectra(h).numColdSpectra = sum(vertcat(calibratedSpectra(goodSpectra).numColdSpectra));
        integratedSpectra(h).numSkySpectra = sum(vertcat(calibratedSpectra(goodSpectra).numAntSpectra));
        
        meanAngleAT = nanmean([calibratedSpectra(goodSpectra).meanAngleAntenna]);
        tod = nanmean([calibratedSpectra(goodSpectra).TOD]);
        dateTime=nanmean([calibratedSpectra(goodSpectra).dateTime]);
    end
    
    integratedSpectra(h).TSys=nanmean([calibratedSpectra(goodSpectra).TSys]);
    if isnan(integratedSpectra(h).TSys)
        integratedSpectra(h).TSys=-9999;
    end
    
    integratedSpectra(h).THot=nanmean([calibratedSpectra(goodSpectra).THot]);
    if isnan(integratedSpectra(h).THot)
        integratedSpectra(h).THot=-9999;
    end
    
    %stdTHot
    %StdTSYS...
    
    integratedSpectra(h).TWindow=nanmean([calibratedSpectra(indSpectra).TWindow]);
    if isfield(calibratedSpectra,'TOut')
        integratedSpectra(h).TOut=nanmean([calibratedSpectra(indSpectra).TOut]);
    else
        integratedSpectra(h).TOut = -9999;
    end

    integratedSpectra(h).if=calibratedSpectra(indSpectra(1)).if;
    integratedSpectra(h).freq=calibratedSpectra(indSpectra(1)).freq;
    integratedSpectra(h).timeMin=calibratedSpectra(indSpectra(1)).timeMin;
    
    integratedSpectra(h).firstSkyTime=calibratedSpectra(indSpectra(1)).firstSkyTime;
    integratedSpectra(h).lastSkyTime=calibratedSpectra(indSpectra(end)).lastSkyTime;
    
    integratedSpectra(h).year=calibratedSpectra(indSpectra(1)).year;
    integratedSpectra(h).month=calibratedSpectra(indSpectra(1)).month;
    integratedSpectra(h).day=calibratedSpectra(indSpectra(1)).day;
    
    integratedSpectra(h).calibrationTime=calibratedSpectra(1).calibrationTime;
    integratedSpectra(h).integrationTime=calibrationTool.integrationTime*60;
    
    % variable that we want to integrate with good spectra if exist
    integratedSpectra(h).Tb=integratedTb;
    integratedSpectra(h).stdTb=integratedStdTb;
    integratedSpectra(h).meanStdTbFromCal=integratedMeanStdTbFromCalMean;
    integratedSpectra(h).meanAngleAntenna = meanAngleAT;
    integratedSpectra(h).TOD=tod;
    integratedSpectra(h).dateTime=dateTime;
    
    % Meteo Data are integrated on all calibrated spectra
    integratedSpectra(h).meanAirPressure=nanmean([calibratedSpectra(indSpectra).meanAirPressure]);
    if isnan(integratedSpectra(h).meanAirPressure)
        integratedSpectra(h).meanAirPressure=-9999;
    end
    
    integratedSpectra(h).meanAirTemperature=nanmean([calibratedSpectra(indSpectra).meanAirTemperature]);   
    if isnan(integratedSpectra(h).meanAirTemperature)
        integratedSpectra(h).meanAirTemperature=-9999;
    end
    
    integratedSpectra(h).meanRelativeHumidity=nanmean([calibratedSpectra(indSpectra).meanRelHumidity]);
    if isnan(integratedSpectra(h).meanRelativeHumidity)
        integratedSpectra(h).meanRelativeHumidity=-9999;
    end
    
    integratedSpectra(h).rainAccumulation=nansum([calibratedSpectra(indSpectra).rainAccumulation]);
    if isnan(integratedSpectra(h).rainAccumulation)
        integratedSpectra(h).rainAccumulation=-9999;
    end
    
end
level1b.integration=integratedSpectra;
end
