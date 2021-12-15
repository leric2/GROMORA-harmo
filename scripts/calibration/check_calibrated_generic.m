function [calibratedSpectra, logFile] = check_calibrated_generic(logFile,calibrationTool,calibratedSpectra)
%==========================================================================
% NAME      | check_calibrated_generic.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2020
%           |
% ABSTRACT  | Function performing some checks on the calibration and
%           | adding meta information to the calibratedSpectra structure.
%           | It builts also the flags vector for the level1a data for each
%           | calibration cycle and add every new information to the
%           | calibrated spectra structure (IN/OUT).
%           | 
%           |
% ARGUMENTS | INPUTS:   1. logFile: harmonized log file 
%           |           2. calibratedSpectra
%           |           3. calibrationTool:
%           |               - minNumberOfIndicePerCycle
%           |               - IQProcessing
%           |               - samplingRateFFTS
%           |               - numberOfChannels
%           |               - observationFreq
%           |               - LOFreqTot
%           |               - frequencyBandAroundCenterTNoise
%           |               - TCold
%           |               - TNoiseCenterTh
%           |               - TNoiseThresh
%           |               - stdTNoise
%           |               - stdTNoiseThresh
%           |               - maxStdDevTbCal
%           |               - maxProportionOfIndLN2SensorOutlier
%           |               - maxProportionOfIndLN2LevelOutlier
%           |               - hotTemperatureStdThreshold
%           |               - maxProportionOfIndFFTadcOverload
%           |               - numberOfAquisitionSpectraAntenna
%           |               - numberOfAquisitionSpectraHot
%           |               - numberOfAquisitionSpectraCold
%           |               - backgroundMWTb
%           |               - TC
%           |               - stdAntAngleThresh
%           |               - referenceTime
%           | 
%           | OUTPUTS:  1. calibratedSpectra
%           |
%==========================================================================

% Checking all calibration cycle
for i = 1:size(calibratedSpectra,2)
    ia=calibratedSpectra(i).antennaIndCleanAngle;
    ih=calibratedSpectra(i).hotInd;
    ic=calibratedSpectra(i).coldInd;
    
    ind=sort([ia ih ic]);
    
    %%%%%%%%%%% Flag 1 %%%%%%%%%%%
    % The number of indices for the 3 positions:
    if ((length(ih) > calibrationTool.minNumberOfIndicePerCycle) &&...
            (length(ia) > calibrationTool.minNumberOfIndicePerCycle) &&...
            (length(ic) > calibrationTool.minNumberOfIndicePerCycle))
        sufficientNumberOfIndices=1;
    else
        sufficientNumberOfIndices=0;
    end
    
    %calibratedSpectra(i).TbPlanck =  (calibrationTool.h*calibratedSpectra(i).freq/calibrationTool.kb)./log((2*calibrationTool.h*calibratedSpectra(i).freq.^3)./(calibratedSpectra(i).intensityPlanck*calibrationTool.lightSpeed^2) + 1);
    % from TbRJE: same results in TbPlanck as from intensity calibration !
    %I_RJE =  calibratedSpectra(i).TbRJE.*(2*calibrationTool.kb*calibratedSpectra(i).freq.^2)/calibrationTool.lightSpeed^2;
    %calibratedSpectra(i).TbPlanck2 = (calibrationTool.h*calibratedSpectra(i).freq/calibrationTool.kb)./log((2*calibrationTool.h*calibratedSpectra(i).freq.^3)./(I_RJE*calibrationTool.lightSpeed^2) + 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Noise Receiver Temperature
    % Computing TN around the line center (approximately +- 200 MHz)
    centerChannels=find(calibratedSpectra(i).freq>=calibratedSpectra(i).observationFreq-calibrationTool.frequencyBandAroundCenterTNoise & calibratedSpectra(i).freq<calibratedSpectra(i).observationFreq+calibrationTool.frequencyBandAroundCenterTNoise);
    
    Ycenter=calibratedSpectra(i).Yspectral(centerChannels);
    
    % Removing extreme outliers before computing TN:
    boxCarFilter=ones(100,1)/100;
    Yfiltered=conv(Ycenter,boxCarFilter,'same');
    
    YcenterSmoothed=Ycenter(abs(Ycenter-Yfiltered)<2*nanstd(Ycenter));
    
    TNoiseCenter=(calibratedSpectra(i).THot-YcenterSmoothed*calibrationTool.TCold)./(YcenterSmoothed-1);
    
    %TSysCenter2=(calibratedSpectra(i).THot-mean(Ycenter)*retrievalTool.TCold)./(mean(Ycenter)-1);
    
    calibratedSpectra(i).TNoise=nanmean(TNoiseCenter);
    
    %%%%%%%%%%% Flag 2 %%%%%%%%%%%
    % stdTNoise defined on drift !!!!???
    if (abs(calibratedSpectra(i).TNoise-calibrationTool.TNoiseCenterTh)>calibrationTool.TNoiseThresh | calibratedSpectra(i).stdTNoise > calibrationTool.stdTNoiseThresh)
        noiseTemperatureOK=0;
    else
        noiseTemperatureOK=1;
    end
    
    % TNoise from the log file
    calibratedSpectra(i).TNoiseLog=nanmean(logFile.FE_T_Sys(ind));
    calibratedSpectra(i).stdTNoiseLog=nanstd(logFile.FE_T_Sys(ind));
    
    % Mean standard deviation of calibrated sky measurement during the
    % cycle
    calibratedSpectra(i).potentialBadChannels = calibratedSpectra(i).stdTb > calibrationTool.maxStdDevTbCal;
    calibratedSpectra(i).meanStdTb=nanmean(calibratedSpectra(i).stdTb(~calibratedSpectra(i).potentialBadChannels));
    
    diff_Tb = diff(calibratedSpectra(i).Tb(~calibratedSpectra(i).potentialBadChannels));
    calibratedSpectra(i).noiseLevel = nanstd(diff_Tb(~isinf(diff_Tb)))/sqrt(2);

    % Mean calibrated Tb for this cycle:
    calibratedSpectra(i).meanTb = nanmean(calibratedSpectra(i).Tb(~calibratedSpectra(i).potentialBadChannels));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cold and Hot spectra variation among this cycle
    calibratedSpectra(i).meanStdHotSpectra=nanmean(calibratedSpectra(i).stdHotSpectra);
    calibratedSpectra(i).meanStdColdSpectra=nanmean(calibratedSpectra(i).stdColdSpectra);
    
    % Mean counts
    calibratedSpectra(i).meanHotCounts = nanmean(calibratedSpectra(i).meanHotSpectra);
    calibratedSpectra(i).meanColdCounts = nanmean(calibratedSpectra(i).meanColdSpectra);
    calibratedSpectra(i).meanSkyCounts = nanmean(calibratedSpectra(i).meanSkySpectra);
    
    %%%%%%%%%%% Flag 3 %%%%%%%%%%%    
    % Liquid Nitrogen sensors
    if sum(~logFile.LN2_Sensors_OK(ind)==1) > calibrationTool.maxProportionOfIndLN2SensorOutlier*length(ind)
        LN2SensorsOK=0;
    else
        LN2SensorsOK=1;
    end
    
    %%%%%%%%%%% Flag 4 %%%%%%%%%%%
    % Liquid Nitrogen level
    if sum(~logFile.LN2_Level_OK(ind)==1) > calibrationTool.maxProportionOfIndLN2LevelOutlier*length(ind)
        LN2LevelOK=0;
    else
        LN2LevelOK=1;
    end

    %%%%%%%%%%% Flag 5 %%%%%%%%%%%
    % Hot load check flag
    % abs(calibratedSpectra(i).THot - calibrationTool.THotTh)>calibrationTool.THotAbsThresh |
    if (calibratedSpectra(i).stdTHot>calibrationTool.hotTemperatureStdThreshold)
        hotLoadOK=0;
    else
        hotLoadOK=1;
    end
    
    %%%%%%%%%% Flag 6 %%%%%%%%%%%
    %FFTS aquisition variable (not used for now but might be useful
    calibrationTool.maxProportionOfIndFFTadcOverload = 0.1;
    if sum(~(logFile.FFT_adc_overload(ind) == 0)) > calibrationTool.maxProportionOfIndFFTadcOverload*length(ind)
        FFT_adc_overload_OK=0;
    else
        FFT_adc_overload_OK=1;
    end
    
   
    %%%%%%%%%% Flag 6 %%%%%%%%%%%
    calibratedSpectra(i).stdVGunn = std(logFile.V_Gunn(ind),'omitnan');
    calibratedSpectra(i).VGunn = mean(logFile.V_Gunn(ind),'omitnan');
    if (calibratedSpectra(i).stdVGunn > calibrationTool.maxStdV_Gun) | (sum((logFile.Freq_Lock(ind) == 0)) > calibrationTool.maxProportionFreqLockError*length(ind)) 
        freq_lock_OK=0;
    else
        freq_lock_OK=1;
    end
    
    
    %%%%%%%%%%% Flag ... %%%%%%%%%%%    NOT USED 
    % FFTS number of aquisition
    if ((all(logFile.FFT_Nr_of_acq(ia))==calibrationTool.numberOfAquisitionSpectraAntenna) & (all(logFile.FFT_Nr_of_acq(ih))==calibrationTool.numberOfAquisitionSpectraHot) & (all(logFile.FFT_Nr_of_acq(ih))==calibrationTool.numberOfAquisitionSpectraCold))
        calibratedSpectra(i).FFT_nr_aquisition_OK=1;
    else
        calibratedSpectra(i).FFT_nr_aquisition_OK=0;
    end
    
    %%%%%%%%%%% Flag 7 %%%%%%%%%%%
    % Antenna angle
    calibratedSpectra(i).meanAngleAntenna=mean(logFile.Elevation_Angle(ia));
    calibratedSpectra(i).stdAngleAntenna=std(logFile.Elevation_Angle(ia));
    
    if calibratedSpectra(i).stdAngleAntenna > calibrationTool.stdAntAngleThresh
        PointingAngleOK=0;
    else
        PointingAngleOK=1;
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding other variables (if existing)
    if isfield(logFile,'T_Room')
        calibratedSpectra(i).TempRoom=nanmean(logFile.T_Room(ind));
        calibratedSpectra(i).stdTempRoom=nanstd(logFile.T_Room(ind));
    else
        calibratedSpectra(i).TempRoom=-9999;
        calibratedSpectra(i).stdTempRoom=-9999;
    end
    if isfield(logFile,'T_Out')
        calibratedSpectra(i).TempOut=nanmean(logFile.T_Out(ind));
        calibratedSpectra(i).stdTempOut=nanstd(logFile.T_Out(ind));
    else
        calibratedSpectra(i).TempOut=-9999;
        calibratedSpectra(i).stdTempOut=-9999;
    end
    if isfield(logFile,'T_Window')
        calibratedSpectra(i).TempWindow=nanmean(logFile.T_Window(ind));
        calibratedSpectra(i).stdTempWindow=nanstd(logFile.T_Window(ind));
    else
        calibratedSpectra(i).TempWindow=-9999;
        calibratedSpectra(i).stdTempWindow=-9999;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time variable for this cycle
    
    % As we are always using daily raw files:
    calibratedSpectra(i).year=calibratedSpectra(i).theoreticalStartTime.Year;
    calibratedSpectra(i).month=calibratedSpectra(i).theoreticalStartTime.Month;
    calibratedSpectra(i).day=calibratedSpectra(i).theoreticalStartTime.Day;
    
    calibratedSpectra(i).timeMin = datestr(calibratedSpectra(i).theoreticalStartTime,'YYYY_mm_dd_HH:MM:SS');
    calibratedSpectra(i).timeMax = datestr(calibratedSpectra(i).theoreticalStartTime + minutes(calibratedSpectra(i).calibrationTime),'YYYY_mm_dd_HH:MM:SS');
    
    calibratedSpectra(i).timeMin=datenum(calibratedSpectra(i).timeMin,'YYYY_mm_dd_HH:MM:SS')-calibrationTool.referenceTime;
    calibratedSpectra(i).timeMax=datenum(calibratedSpectra(i).timeMax,'YYYY_mm_dd_HH:MM:SS')-calibrationTool.referenceTime;
 
    if ~isempty(ia)
        %calibratedSpectra(i).dateStart=datestr(logFile.x(1:6,ia(1))','yyyymmddTHHMMSSZ');
        %calibratedSpectra(i).dateStop=datestr(logFile.x(1:6,ia(end))','yyyymmddTHHMMSSZ');
        
        %calibratedSpectra(i).firstSkyTime=datenum(calibratedSpectra(i).dateStart,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
        %calibratedSpectra(i).lastSkyTime=datenum(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
        
        %only possible if ia is not empty !
        calibratedSpectra(i).firstSkyTime=logFile.time(ia(1))-calibrationTool.referenceTime;
        calibratedSpectra(i).lastSkyTime=logFile.time(ia(end))-calibrationTool.referenceTime;
        
        % "mean time" of the calibration cycle (mean of all antenna measurements)
        calibratedSpectra(i).meanAntTime = nanmean(logFile.dateTime(ia));
        calibratedSpectra(i).meanAntTimeUTC = nanmean(logFile.dateTime(ia));
        calibratedSpectra(i).meanAntTimeUTC.TimeZone='Z';
        
        calibratedSpectra(i).meanDatetime = datenum(calibratedSpectra(i).meanAntTime)-calibrationTool.referenceTime;
        
        calibratedSpectra(i).meanDatetimeMJD2K = convert_to_MJD2K(...
            calibratedSpectra(i).meanAntTimeUTC.Year,...
            calibratedSpectra(i).meanAntTimeUTC.Month,...
            calibratedSpectra(i).meanAntTimeUTC.Day,...
            calibratedSpectra(i).meanAntTimeUTC.Hour,...
            calibratedSpectra(i).meanAntTimeUTC.Minute,...
            calibratedSpectra(i).meanAntTimeUTC.Second);
        
        %calibratedSpectra(i).timeOfDay = 24*(datenum(calibratedSpectra(i).meanAntTime) -datenum(calibratedSpectra(i).year,calibratedSpectra(i).month,calibratedSpectra(i).day));
        
        calibratedSpectra(i).timeOfDay = calibratedSpectra(i).meanAntTime.Hour + (calibratedSpectra(i).meanAntTime.Minute*60 + calibratedSpectra(i).meanAntTime.Second)/3600;
        %calibratedSpectra(i).meanDatetimeUnit='days since 1970-01-01 00:00:00';
        %calibratedSpectra(i).calendar='standard';
    
        %calibratedSpectra(i).startTimeInt8=int8(datestr(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ'));
    
    else
        calibratedSpectra(i).firstSkyTime = -9999;
        calibratedSpectra(i).lastSkyTime = -9999;
        calibratedSpectra(i).timeOfDay = -9999;
        
        
        calibratedSpectra(i).meanAntTime = calibratedSpectra(i).theoreticalStartTime + 0.5*minutes(calibratedSpectra(i).calibrationTime);
        calibratedSpectra(i).meanAntTimeUTC = calibratedSpectra(i).theoreticalStartTime + 0.5*minutes(calibratedSpectra(i).calibrationTime);
        calibratedSpectra(i).meanAntTimeUTC.TimeZone='Z';
        calibratedSpectra(i).meanDatetime = datenum(calibratedSpectra(i).meanAntTime)-calibrationTool.referenceTime;
        
        calibratedSpectra(i).meanDatetimeMJD2K = convert_to_MJD2K(...
            calibratedSpectra(i).meanAntTimeUTC.Year,...
            calibratedSpectra(i).meanAntTimeUTC.Month,...
            calibratedSpectra(i).meanAntTimeUTC.Day,...
            calibratedSpectra(i).meanAntTimeUTC.Hour,...
            calibratedSpectra(i).meanAntTimeUTC.Minute,...
            calibratedSpectra(i).meanAntTimeUTC.Second);
    end
    
    if isfield(logFile,'TC')
       isTC = isbetween([logFile.TC.skyMeanDatetime],calibratedSpectra(i).theoreticalStartTime,calibratedSpectra(i).theoreticalStartTime+minutes(calibrationTool.calibrationTime));
       % isTC = ([logFile.TC.meanDateNum] >= calibratedSpectra(i).firstSkyTime & [logFile.TC.meanDateNum] < calibratedSpectra(i).lastSkyTime);
       % check if there was a tc done during this cycle (only mean
       % datetime)
       if sum(isTC) > 0
           logFile.TC(isTC).coldSpectra = calibratedSpectra(i).meanColdSpectra(logFile.TC(isTC).channels);
           logFile.TC(isTC).coldCalib = nanmean(logFile.TC(isTC).coldSpectra);
           logFile.TC(isTC).hotSpectra = calibratedSpectra(i).meanHotSpectra(logFile.TC(isTC).channels);
           logFile.TC(isTC).hotCalib = nanmean(logFile.TC(isTC).hotSpectra);
           logFile.TC(isTC).THotCalib = calibratedSpectra(i).THot;
           logFile.TC(isTC).meanFreq = mean(calibratedSpectra(i).freq(logFile.TC(isTC).channels));
           logFile.TC(isTC).frequency = calibratedSpectra(i).freq(logFile.TC(isTC).channels);
           
           if isfield(logFile.meteo,'air_temperature')
               air_temp = [logFile.meteo.air_temperature];
               % we take the 3 meteo data aroung this calibration cycle to have
               % at least one value.
               meteoInd = isbetween([logFile.meteo.dateTime], calibratedSpectra(i).theoreticalStartTime-minutes(calibrationTool.calibrationTime),calibratedSpectra(i).theoreticalStartTime+minutes(calibrationTool.calibrationTime));
               % just for estimation
               Teff = mean(air_temp(meteoInd))-calibrationTool.TC.deltaT;
           else
               disp('we said, no meteo data found so lets make a guess for Tair (10 degC)');
               Teff = 283 - calibrationTool.TC.deltaT; 
           end
           logFile.TC(isTC).Tb_Calib = calibrationTool.TCold + (logFile.TC(isTC).THotCalib - calibrationTool.TCold) .* (logFile.TC(isTC).sky - logFile.TC(isTC).coldCalib)./(logFile.TC(isTC).hotCalib - logFile.TC(isTC).coldCalib);
           
           % In terms of RJE:
           TbColdRJ = planck_function(calibrationTool, calibrationTool.TCold, logFile.TC(isTC).frequency)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*( logFile.TC(isTC).frequency.^2));
           TbHotRJ = planck_function(calibrationTool, logFile.TC(isTC).THotCalib, logFile.TC(isTC).frequency)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*( logFile.TC(isTC).frequency.^2));
           logFile.TC(isTC).TbRJE = TbColdRJ + (TbHotRJ - TbColdRJ) .* (logFile.TC(isTC).sky_spectra - logFile.TC(isTC).coldSpectra)./(logFile.TC(isTC).hotSpectra - logFile.TC(isTC).coldSpectra);
                
           TC(isTC).TbRJE_mean = nanmean(logFile.TC(isTC).TbRJE,2);

           TmeanRJ = planck_function(calibrationTool, Teff, logFile.TC(isTC).meanFreq)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*(logFile.TC(isTC).meanFreq^2));
           TbgRJ = planck_function(calibrationTool, calibrationTool.backgroundMWTb, logFile.TC(isTC).meanFreq)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*(logFile.TC(isTC).meanFreq^2)); 
                
           tau_slantRJ = log((TmeanRJ-TbgRJ)./(TmeanRJ-TC(isTC).TbRJE_mean));                
           
           tau_slant = log((Teff-calibrationTool.backgroundMWTb)./(Teff-logFile.TC(isTC).Tb_Calib));
           am = 1./sind(logFile.TC(isTC).skyAngle);
           
           % fit the airmass-slant opacity data pairs (RJE)
           [pRJE,s] = polyfit (am, tau_slantRJ, 1);
           logFile.TC(isTC).tauCalibZenith = pRJE(1);
           
           % fit the airmass-slant opacity data pairs
           %[p,s] = polyfit (am, tau_slant, 1);
           %logFile.TC(isTC).tauCalibZenith = p(1);
           logFile.TC(isTC).tauCalib = logFile.TC(isTC).tauCalibZenith * 1/sind(calibratedSpectra(i).meanAngleAntenna);
       end
       
    end
%     % Effective calibration time for this cycle (TO CHECK IF NEEDED ?)
%     effectiveTime=0;
%     %if i<size(calibratedSpectra,2)
%     for l= 1:length(ind)-1
%         effectiveTime=effectiveTime+(logFile.t(ind(l+1))-logFile.t(ind(l)))*60;
%     end
%     % adding the last time if possible
%     if i<size(calibratedSpectra,2)
%         effectiveTime=effectiveTime+(logFile.t(ind(l+1)+1)-logFile.t(ind(l+1)));
%     else
%         % TODO
%         effectiveTime=effectiveTime;
%     end
   % calibratedSpectra(i).effectiveCalibrationTime=effectiveTime;
      
    %calibratedSpectra(i).effectiveCalibrationTimeHot=length(ih)*retrievalTool.calibTimeHot;
    %calibratedSpectra(i).effectiveCalibrationTimeAntenna=length(ia)*retrievalTool.calibTimeAntenna;
    %calibratedSpectra(i).effectiveCalibrationTimeCold=length(ic)*retrievalTool.calibTimeCold;
    %calibratedSpectra(i).effectiveCalibrationTime=calibratedSpectra(i).effectiveCalibrationTimeHot+calibratedSpectra(i).effectiveCalibrationTimeAntenna+calibratedSpectra(i).effectiveCalibrationTimeCold;
   
    calibratedSpectra(i).numberOfIndices=[
        length(ih),...
        length(ic),...
        length(ia)];
    
    % Error vector for this calibration cycle
    calibratedSpectra(i).errorVector=[
        sufficientNumberOfIndices,...
        noiseTemperatureOK,...
        LN2SensorsOK,...
        LN2LevelOK,...
        hotLoadOK,...
        PointingAngleOK,...
        freq_lock_OK];
    
    % Error vector description:
    calibratedSpectra(i).errorVectorDescription=[
        "sufficientNumberOfIndices",...
        "noiseTemperatureOK",...
        "LN2SensorsOK",...
        "LN2LevelOK",...
        "hotLoadOK",...
        "pointingAngleOK",...
        "FreqLockOK"];


    calibratedSpectra(i).outlierCalib = NaN;
    warning('off','backtrace')
    if (sum(calibratedSpectra(i).errorVector)<length(calibratedSpectra(i).errorVector))
        calibratedSpectra(i).outlierCalib = 1;
        errorV=num2str(calibratedSpectra(i).errorVector);
        %disp(['Calibration Cycle number ' num2str(i) ', TOD: ' datestr(timeofday(calibratedSpectra(i).meanAntTime),'HH:MM:SS')])
        disp(['Problem with calibration n. ' num2str(i) ', TOD: ' datestr(timeofday(calibratedSpectra(i).meanAntTime),'HH:MM:SS') ', error: ']);
        disp(errorV)
        disp(calibratedSpectra(i).errorVectorDescription(~calibratedSpectra(i).errorVector))
    end
%     if strcmp(calibrationTool.outlierDectectionType,'RFI')
%         if calibratedSpectra(i).outlierRFISky + calibratedSpectra(i).outlierDetectColdRFI  + calibratedSpectra(i).outlierDetectHotRFI  > 0
%             disp(['Potential RFI problem n. ' num2str(i) ', TOD: ' datestr(timeofday(calibratedSpectra(i).meanAntTime),'HH:MM:SS')]);
%             fig=figure("visible","off");
%             ax1 = subplot(3,1,1); plot(ax1, calibratedSpectra(i).freqRFI,calibratedSpectra(i).AntSpectraRFI); ylim([-20, 50]); title('Sky')
%             ax2 = subplot(3,1,2); plot(ax2, calibratedSpectra(i).freqRFI,calibratedSpectra(i).ColdSpectraRFI);ylim([-20, 50]); title('Cold')
%             ax3 = subplot(3,1,3); plot(ax3, calibratedSpectra(i).freqRFI,calibratedSpectra(i).HotSpectraRFI);ylim([-20, 50]); title('Hot'); xlabel('IF [MHz]')
%             
%             print(fig,['/home/esauvageat/Desktop/RFI/GROMOS_RFI_' calibrationTool.dateStr '_' num2str(i)],'-dpdf','-fillpage')
%             xlabel('IF [MHz]')
%             ylabel('count difference')
%         end
%     end
   
end

end

