function integratedSpectra = integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra)
%==========================================================================
% NAME          | integrate_calibrated_spectra_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 03.2020
%               |
% ABSTRACT      | Function for integrating calibrated spectra from an
%               | instrument. It takes the standard calibratedSpectra
%               | structure and performs the integration on the defined 
%               | integrationTime. Note that for some variables (e.g. for 
%               | precipitation), this is more than a simple averaging. 
%               | The results is saved into integratedSpectra, a similar
%               | structure to calibratedSpectra but on a timescale 
%               | corresponding to the integration times. 
%               |
% ARGUMENTS     | INPUTS:  1. calibrationTool:
%               |           - integrationTime
%               |           - Year, Month, Day
%               |           - timeZone
%               |           - integrationTime
%               |           - filterByTransmittance
%               |           - transmittanceThreshold
%               |           - filterByFlags
%               |           - flagVectorLength
%               |           - numberOfChannels
%               |          2. calibratedSpectra
%               |
%               | OUTPUTS: 1. integratedSpectra
%               |
%==========================================================================
% Threshold for the integration, calibTime has to be in [min]
dt=hours(calibrationTool.integrationTime/60);
%timeThresh=0:dt:23;
timeThresh = datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day):dt:datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day+1);
timeThresh.TimeZone = calibrationTool.timeZone;

for h = 1:length(timeThresh)-1
    % Finding the spectra during this time stamp:
    indSpectra=find([calibratedSpectra.dateTime]>=timeThresh(h) & [calibratedSpectra.dateTime]<timeThresh(h)+dt);
    
    % Selecting only spectra with:
    
    % transmission > 0.2 in the integration:
    if calibrationTool.filterByTransmittance & ~calibrationTool.missing_meteo
        goodSpectra=indSpectra([calibratedSpectra(indSpectra).troposphericTransmittance] >= calibrationTool.transmittanceThreshold);
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
        integratedPlanckIntensity=-9999*ones(1,calibrationTool.numberOfChannels);
        integratedStdTb = -9999*ones(1,calibrationTool.numberOfChannels);
        integratedMeanStdTbFromCalMean = -9999;
        integratedSpectra(h).numberOfAveragedSpectra = 0;
        % goodSpectra=indSpectra;
        
        integratedSpectra(h).number_of_hot_spectra = 0;
        integratedSpectra(h).number_of_cold_spectra = 0;
        integratedSpectra(h).number_of_sky_spectra = 0;
        
        integratedSpectra(h).tropospheric_opacity_tc = -9999;
        meanAngleAT = -9999;
        tod = nanmean([calibratedSpectra(indSpectra).time_of_day]);
        dateTime=nanmean([calibratedSpectra(indSpectra).dateTime]);
    else
        integratedSpectra(h).numberOfAveragedSpectra=length(goodSpectra);
        % Averaging the good spectra together
        integratedTb=nanmean(vertcat(calibratedSpectra(goodSpectra).Tb),1);
        if calibrationTool.savePlanckIntensity
            integratedPlanckIntensity = nanmean(vertcat(calibratedSpectra(goodSpectra).intensity_planck),1);
        end
        % Computing the std deviation on the integrated spectra
        sum_of_variance = nansum(vertcat(calibratedSpectra(goodSpectra).stdTb).^2,1);
        %integratedStdTb = nansum(vertcat(calibratedSpectra(goodSpectra).stdTb),1) / sqrt(length(goodSpectra));
        integratedStdTb = sqrt(sum_of_variance/length(goodSpectra));
        %integratedMeanStdTbFromCalMean = sqrt(nansum([calibratedSpectra(goodSpectra).meanStdTb].^2)/length(goodSpectra))
        
        % For the mean, we take only good quality channels
        good_quality_all = any(vertcat(calibratedSpectra(goodSpectra).channelsQuality),1);
        integratedMeanStdTbFromCalMean = mean(integratedStdTb(good_quality_all));
        
        % Summing the number of spectra for hot, cold and antenna:
        integratedSpectra(h).number_of_hot_spectra = sum(vertcat(calibratedSpectra(goodSpectra).number_of_hot_spectra));
        integratedSpectra(h).number_of_cold_spectra = sum(vertcat(calibratedSpectra(goodSpectra).number_of_cold_spectra));
        integratedSpectra(h).number_of_sky_spectra = sum(vertcat(calibratedSpectra(goodSpectra).number_of_sky_spectra));
        
        integratedSpectra(h).tropospheric_opacity_tc = nanmean([calibratedSpectra(goodSpectra).troposphericOpacityTC]);
        meanAngleAT = nanmean([calibratedSpectra(goodSpectra).mean_sky_elevation_angle]);
        tod = nanmean([calibratedSpectra(goodSpectra).time_of_day]);
        dateTime=nanmean([calibratedSpectra(goodSpectra).dateTime]);
    end
    
    integratedSpectra(h).TNoise=nanmean([calibratedSpectra(goodSpectra).TNoise]);
    if isnan(integratedSpectra(h).TNoise)
        integratedSpectra(h).TNoise=-9999;
    end
    
    integratedSpectra(h).THot=nanmean([calibratedSpectra(goodSpectra).THot]);
    if isnan(integratedSpectra(h).THot)
        integratedSpectra(h).THot=-9999;
    end
    
    if isnan(integratedSpectra(h).tropospheric_opacity_tc)
        integratedSpectra(h).tropospheric_opacity_tc = -9999;
    end
    
    %stdTHot
    %StdTNoise...
    
    integratedSpectra(h).TWindow=nanmean([calibratedSpectra(indSpectra).TWindow]);
    integratedSpectra(h).TRoom=nanmean([calibratedSpectra(indSpectra).TRoom]);
    integratedSpectra(h).VGunn=nanmean([calibratedSpectra(indSpectra).VGunn]);

    if isfield(calibratedSpectra,'TOut')
        integratedSpectra(h).TOut=nanmean([calibratedSpectra(indSpectra).TOut]);
    else
        integratedSpectra(h).TOut = -9999;
    end

    integratedSpectra(h).intermediate_freq=calibratedSpectra(indSpectra(1)).intermediate_freq;
    integratedSpectra(h).frequencies=calibratedSpectra(indSpectra(1)).frequencies;
    integratedSpectra(h).time_min=calibratedSpectra(indSpectra(1)).time_min;
    
    integratedSpectra(h).first_sky_time=calibratedSpectra(indSpectra(1)).first_sky_time;
    integratedSpectra(h).last_sky_time=calibratedSpectra(indSpectra(end)).last_sky_time;
    
    integratedSpectra(h).year=calibratedSpectra(indSpectra(1)).year;
    integratedSpectra(h).month=calibratedSpectra(indSpectra(1)).month;
    integratedSpectra(h).day=calibratedSpectra(indSpectra(1)).day;
    
    integratedSpectra(h).calibration_time=calibratedSpectra(1).calibration_time;
    integratedSpectra(h).integration_time=calibrationTool.integrationTime*60;
    
    % variable that we want to integrate with good spectra if exist
    integratedSpectra(h).Tb=integratedTb;
    if calibrationTool.savePlanckIntensity
        integratedSpectra(h).intensity_planck=integratedPlanckIntensity;
    end
    integratedSpectra(h).stdTb=integratedStdTb;
    integratedSpectra(h).meanStdTbFromCal=integratedMeanStdTbFromCalMean;
    integratedSpectra(h).mean_sky_elevation_angle = meanAngleAT;
    integratedSpectra(h).time_of_day=tod;
    integratedSpectra(h).dateTime=dateTime;
    
    % Meteo Data are integrated on all calibrated spectra
    integratedSpectra(h).mean_air_pressure=nanmean([calibratedSpectra(indSpectra).mean_air_pressure]);
    if isnan(integratedSpectra(h).mean_air_pressure)
        integratedSpectra(h).mean_air_pressure=-9999;
    end
    
    integratedSpectra(h).mean_air_temperature=nanmean([calibratedSpectra(indSpectra).mean_air_temperature]);   
    if isnan(integratedSpectra(h).mean_air_temperature)
        integratedSpectra(h).mean_air_temperature=-9999;
    end
    
    integratedSpectra(h).mean_relative_humidity=nanmean([calibratedSpectra(indSpectra).mean_relative_humidity]);
    if isnan(integratedSpectra(h).mean_relative_humidity)
        integratedSpectra(h).mean_relative_humidity=-9999;
    end
    
    integratedSpectra(h).rain_accumulation=nansum([calibratedSpectra(indSpectra).rain_accumulation]);
    if isnan(integratedSpectra(h).rain_accumulation)
        integratedSpectra(h).rain_accumulation=-9999;
    end
end
end
