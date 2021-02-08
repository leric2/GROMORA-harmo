function integratedSpectra = check_integrated_generic(calibrationTool,integratedSpectra)
%==========================================================================
% NAME      | check_integrated_generic.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2020
%           |
% ABSTRACT  | Function performing some checks on the integration and
%           | adding meta information to the integratedSpectra structure.
%           | It builts also the flags vector for the level1b data for each
%           | integration cycle and add every new information to the
%           | integrated spectra structure (IN/OUT).
%           | 
%           |
% ARGUMENTS | INPUTS:   1. integratedSpectra: standard structure containing
%           |             the integrated data.
%           |           2. calibrationTool:
%           |               - referenceTime
%           |               - integrationTime
%           |               - cycleDurationCold
%           |               - cycleDurationHot
%           |               - cycleDurationSky
%           |               - troposphericCorrection
%           |               - minNumberOfAvgSpectra
%           |               - rainAccumulationThreshold
%           |               - troposphericTransmittanceThresh
%           |
%           | OUTPUTS:  1. integratedSpectra
%           |
%==========================================================================
% Checking all integration cycles
for i = 1:size(integratedSpectra,2)
    
    integratedSpectra(i).meanDatetime=datenum(integratedSpectra(i).dateTime)-calibrationTool.referenceTime;
    integratedSpectra(i).meanDatetimeUTC = integratedSpectra(i).dateTime;
    integratedSpectra(i).meanDatetimeUTC.TimeZone='Z';
    
    integratedSpectra(i).MJD2K = convert_to_MJD2K(...
        integratedSpectra(i).meanDatetimeUTC.Year,...
        integratedSpectra(i).meanDatetimeUTC.Month,...
        integratedSpectra(i).meanDatetimeUTC.Day,...
        integratedSpectra(i).meanDatetimeUTC.Hour,...
        integratedSpectra(i).meanDatetimeUTC.Minute,...
        integratedSpectra(i).meanDatetimeUTC.Second);
    
    integratedSpectra(i).numberOfIndices=[
        integratedSpectra(i).number_of_hot_spectra,...
        integratedSpectra(i).number_of_cold_spectra,...
        integratedSpectra(i).number_of_sky_spectra];
    
    if ~isempty(integratedSpectra(i).number_of_cold_spectra)
        integratedSpectra(i).estimatedIntegrationTimeCold = calibrationTool.cycleDurationCold * integratedSpectra(i).number_of_cold_spectra;
    else
        integratedSpectra(i).estimatedIntegrationTimeCold = 0;
    end
    if ~isempty(integratedSpectra(i).number_of_hot_spectra)
        integratedSpectra(i).estimatedIntegrationTimeHot = calibrationTool.cycleDurationHot * integratedSpectra(i).number_of_hot_spectra;
    else
        integratedSpectra(i).estimatedIntegrationTimeHot = 0;
    end
    if ~isempty(integratedSpectra(i).number_of_sky_spectra)
        integratedSpectra(i).estimatedIntegrationTimeSky = calibrationTool.cycleDurationSky * integratedSpectra(i).number_of_sky_spectra;
    else
        integratedSpectra(i).estimatedIntegrationTimeSky = 0;
    end
    if calibrationTool.savePlanckIntensity
        integratedSpectra(i).TbPlanck =  (calibrationTool.h*integratedSpectra(i).frequencies/calibrationTool.kb)./log((2*calibrationTool.h*integratedSpectra(i).frequencies.^3)./(integratedSpectra(i).intensity_planck*calibrationTool.lightSpeed^2) + 1);
        integratedSpectra(i).TbPlanckWinCorr =  (calibrationTool.h*integratedSpectra(i).frequencies/calibrationTool.kb)./log((2*calibrationTool.h*integratedSpectra(i).frequencies.^3)./(integratedSpectra(i).intensityPlanckWinCorr*calibrationTool.lightSpeed^2) + 1);
        integratedSpectra(i).TbPlanckWinTropCorr =  (calibrationTool.h*integratedSpectra(i).frequencies/calibrationTool.kb)./log((2*calibrationTool.h*integratedSpectra(i).frequencies.^3)./(integratedSpectra(i).intensityPlanckTropWindowCorr*calibrationTool.lightSpeed^2) + 1);
    end
    
    %integratedSpectra(i).potentialBadChannels = integratedSpectra(i).stdTb > calibrationTool.maxStdDevTb;
    integratedSpectra(i).potentialBadChannels = isnan(integratedSpectra(i).channelsQuality);
    
    cleanTb = integratedSpectra(i).Tb;
    cleanTb(integratedSpectra(i).potentialBadChannels)=NaN;
    
    if isempty(cleanTb)
        integratedSpectra(i).meanStdTbFromGoodCh = -9999;
        integratedSpectra(i).meanStdFromWings = -9999;
        integratedSpectra(i).meanStdTb = -9999;
        integratedSpectra(i).meanTb =-9999;
    else
        integratedSpectra(i).meanTb = nanmean(cleanTb);
        integratedSpectra(i).meanStdTbFromGoodCh=nanmean(integratedSpectra(i).stdTb(~integratedSpectra(i).potentialBadChannels));
    
        N = length(cleanTb);
        
        skipChannels = calibrationTool.troposphericCorrection.skipFraction * N;
        
        lower = int16(skipChannels);
        upper = int16(skipChannels) + calibrationTool.troposphericCorrection.numberOfChannelsTropCorr;
    
        leftWing = cleanTb(lower:upper);
        rightWing = cleanTb(N-upper:N-lower);
    
        integratedSpectra(i).meanStdFromWings = prctile([leftWing,rightWing],95)-prctile([leftWing,rightWing],5);
        
        diff_Tb = diff(cleanTb);
        integratedSpectra(i).noiseLevel = nanstd(diff_Tb(~isinf(diff_Tb)))/sqrt(2);
    end
    
    %%%%%%%%%%% Flag 1 %%%%%%%%%%%
    % The number of indices for the 3 positions:
    if (integratedSpectra(i).numberOfAveragedSpectra >= calibrationTool.minNumberOfAvgSpectra)
        sufficientNumberOfAvgSpectra=1;
    else
        sufficientNumberOfAvgSpectra=0;
    end
    
    %%%%%%%%%%% Flag 1 %%%%%%%%%%%
    % The number of indices for the 3 positions:
%     if (integratedSpectra(i).numberOfAveragedSpectra > calibrationTool.minNumberOfIndicePerIntegration)
%         sufficientNumberOfIndices=1;
%     else
%         sufficientNumberOfIndices=0;
%         warning('Low number of avg spectra for this integration');
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Antenna angle
    %integratedSpectra(i).stdAngleAntenna=std(standardLog.Elevation_Angle(ia));
    
    %%%%%%%%%%% flag  %%%%%%%%%%
    % Transmittance
    if (integratedSpectra(i).troposphericTransmittance >= calibrationTool.troposphericTransmittanceFlag)
        tropospheric_transmittance_OK=1;
    else
        tropospheric_transmittance_OK=0;
    end
    
    % Rain
    calibrationTool.rainAccumulationThreshold = 0.1*calibrationTool.integrationTime/60;
    if (integratedSpectra(i).rain_accumulation <= calibrationTool.rainAccumulationThreshold)
        rain_Accumulation_OK=1;
    else
        rain_Accumulation_OK=0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error vector for this calibration cycle
    integratedSpectra(i).errorVector=[
        sufficientNumberOfAvgSpectra,...
        tropospheric_transmittance_OK];
    
    % Error vector description:
    integratedSpectra(i).errorVectorDescription=[
        "sufficientNumberOfAvgSpectra",...
        "tropospheric_transmittance_OK"];
    
    integratedSpectra(i).outlierCalib = NaN;

    if (sum(integratedSpectra(i).errorVector)<length(integratedSpectra(i).errorVector))
        integratedSpectra(i).outlierCalib = 1;
        errorV=num2str(integratedSpectra(i).errorVector);
        disp(['Integration n. ' num2str(i) ', TOD: ' datestr(timeofday(integratedSpectra(i).dateTime),'HH:MM:SS')])
        disp(['Problem with this integration, error code: ' errorV]);
        disp(integratedSpectra(i).errorVectorDescription(~integratedSpectra(i).errorVector))
    end
   
end

end

