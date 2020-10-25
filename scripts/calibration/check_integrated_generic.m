function integratedSpectra = check_integrated_generic(calibrationTool,integratedSpectra)
%==========================================================================
% NAME      | check_integrated_generic.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2020
%           |
% ABSTRACT  | Function performing some checks on the calibration and
%           | adding meta information to the integratedSpectra structure.
%           | It builts also the flags vector for the level1a data for each
%           | calibration cycle and add every new information to the
%           | calibrated spectra structure (IN/OUT).
%           | 
%           |
% ARGUMENTS | INPUTS:   - standardLog: harmonized GROSOM log file 
%           |           - calibrationTool
%           |           - integratedSpectra
%           |
%           |
%           | OUTPUTS: - integratedSpectra
%           |
%           |
% CALLS     | 
%           | 
%           | 
%           | 
%           |
%==========================================================================
% Checking all calibration cycle
for i = 1:size(integratedSpectra,2)
    
    integratedSpectra(i).meanDatetime=datenum(integratedSpectra(i).dateTime)-datenum(1970,1,1);
    
    
    
    integratedSpectra(i).numberOfIndices=[
        integratedSpectra(i).numHotSpectra,...
        integratedSpectra(i).numColdSpectra,...
        integratedSpectra(i).numSkySpectra];
    
    integratedSpectra(i).estimatedIntegrationTimeCold = calibrationTool.cycleDurationCold * integratedSpectra(i).numColdSpectra;
    integratedSpectra(i).estimatedIntegrationTimeHot = calibrationTool.cycleDurationHot * integratedSpectra(i).numHotSpectra;
    integratedSpectra(i).estimatedIntegrationTimeSky = calibrationTool.cycleDurationSky * integratedSpectra(i).numSkySpectra;
    
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
    if (integratedSpectra(i).troposphericTransmittance >= 0.2)
        tropospheric_transmittance_OK=1;
    else
        tropospheric_transmittance_OK=0;
    end
    
    % Rain
    calibrationTool.rainAccumulationThreshold = 0.1*calibrationTool.integrationTime/60;
    if (integratedSpectra(i).rainAccumulation <= calibrationTool.rainAccumulationThreshold)
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
        disp(['Calibration Cycle number ' num2str(i) ', TOD: ' datestr(timeofday(integratedSpectra(i).dateTime),'HH:MM:SS')])
        warning(['Problem with this calibration, error code : ' errorV]);
        disp(integratedSpectra(i).errorVectorDescription(~integratedSpectra(i).errorVector))
    end
   
end

end

