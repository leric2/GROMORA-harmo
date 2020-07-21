function calibratedSpectra = check_calibrated_mopi5(logFile,calibrationTool,calibratedSpectra)
%==========================================================================
% NAME          | CHECK_CALIBRATED_SPECTRA Quality check of the calibrated spectra
% TYPE          |
% AUTHOR(S)     |
% CREATION      |
%               |
% ABSTRACT      | Completing and quality checking each calibration cycle
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
% 
%emptyTimestamp=[];
for i =1:size(calibratedSpectra,2)
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
        %warning('Low number of spectra for this cycle');
    end

%     % Effective calibration time for this cycle (TO CHECK IF NEEDED ?)
    %effectiveTime=0;
    %if i<size(calibratedSpectra,2)
%     for l= 1:length(ind)-1
%         effectiveTime=effectiveTime+(log.t(ind(l+1))-log.t(ind(l)))*60;
%     end
%     adding the last time if possible
%     if i<size(calibratedSpectra,2)
%         effectiveTime=effectiveTime+(log.t(ind(l+1)+1)-log.t(ind(l+1)));
%     else
%         TODO
%         effectiveTime=effectiveTime;
%     end
    calibratedSpectra(i).effectiveCalibrationTime=-9999;
      
    %calibratedSpectra(i).effectiveCalibrationTimeHot=length(ih)*retrievalTool.calibTimeHot;
    %calibratedSpectra(i).effectiveCalibrationTimeAntenna=length(ia)*retrievalTool.calibTimeAntenna;
    %calibratedSpectra(i).effectiveCalibrationTimeCold=length(ic)*retrievalTool.calibTimeCold;
    %calibratedSpectra(i).effectiveCalibrationTime=calibratedSpectra(i).effectiveCalibrationTimeHot+calibratedSpectra(i).effectiveCalibrationTimeAntenna+calibratedSpectra(i).effectiveCalibrationTimeCold;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency vector
    if calibrationTool.IQProcessing
        calibratedSpectra(i).if = calibrationTool.samplingRateFFTS/2 * [-1:2/calibrationTool.numberOfChannels:1-2/calibrationTool.numberOfChannels];
        
        calibratedSpectra(i).observationFreq=calibrationTool.observationFreq;
        
        calibratedSpectra(i).LOFreqTot=calibrationTool.LOFreqTot;
    
        calibratedSpectra(i).freq=calibratedSpectra(i).if*1e6+calibratedSpectra(i).LOFreqTot;
    
        calibratedSpectra(i).df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    else
        calibratedSpectra(i).if = calibrationTool.samplingRateFFTS/2 * [0:1/calibrationTool.numberOfChannels:1-1/calibrationTool.numberOfChannels];
    
        calibratedSpectra(i).observationFreq=calibrationTool.observationFreq;
        
        calibratedSpectra(i).LOFreqTot=calibrationTool.LOFreqTot;
    
        calibratedSpectra(i).freq=calibratedSpectra(i).if*1e6+calibratedSpectra(i).LOFreqTot;
    
        calibratedSpectra(i).df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    end
    
    %bw=calibrationTool.instrumentBandwidth;
    %nChannel=calibrationTool.numberOfChannels;
    %df=bw/(nChannel+1); % TOCHECK
    %lc=log.Spectr_line_center(1);

    %calibratedSpectra(i).freq=horzcat(sort(calibratedSpectra(i).LOFreqTot-df*(0:calibrationTool.DCChannel-1)),calibratedSpectra(i).LOFreqTot+df*(1:nChannel-calibrationTool.DCChannel));
    %calibratedSpectra(i).freq=logFile.f;
    
    %calibratedSpectra(i).intFr=calibratedSpectra(i).freq-calibratedSpectra(i).freq(1);
    %calibratedSpectra(i).freq=log.f;
    %calibratedSpectra(i).freq=(calibratedSpectra(i).f0-(lc*df)):df:calibratedSpectra(i).f0+((nChannel-(lc+1))*df);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hot load check flag
    if (calibratedSpectra(i).stdTHot>calibrationTool.hotTemperatureStdThreshold)
        calibratedSpectra(i).hotLoadOK=0;
    else
        calibratedSpectra(i).hotLoadOK=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % System Temperature
    % Computing TN around the line center (approximately +- x MHz)
    
    centerChannels=find(calibratedSpectra(i).freq>=calibratedSpectra(i).observationFreq-calibrationTool.frequencyBandAroundCenterTSys & calibratedSpectra(i).freq<calibratedSpectra(i).observationFreq+calibrationTool.frequencyBandAroundCenterTSys);
    
    Ycenter=calibratedSpectra(i).Yspectral(centerChannels);
    
    % Removing extreme outliers before computing TN:
    boxCarFilter=ones(100,1)/100;
    Yfiltered=conv(Ycenter,boxCarFilter,'same');
    
    YcenterSmoothed=Ycenter(abs(Ycenter-Yfiltered)<2*nanstd(Ycenter));
    
    TSysCenter=(calibratedSpectra(i).THot-YcenterSmoothed*calibrationTool.TCold)./(YcenterSmoothed-1);
    
    %TSysCenter2=(calibratedSpectra(i).THot-mean(Ycenter)*retrievalTool.TCold)./(mean(Ycenter)-1);
    
    calibratedSpectra(i).TSys=nanmean(TSysCenter);
    
    %%%%%%%%%%% Flag 2 %%%%%%%%%%%      
    if (abs(calibratedSpectra(i).TSys-calibrationTool.TSysCenterTh)>calibrationTool.TSysThresh | calibratedSpectra(i).stdTSys > calibrationTool.stdTSysThresh)
        systemTemperatureOK=0;
    else
        systemTemperatureOK=1;
    end
    
    % Tsys from the log file
    calibratedSpectra(i).TSysLog=nanmean(logFile.FE_T_Sys(ind));
    calibratedSpectra(i).stdTSysLog=nanstd(logFile.FE_T_Sys(ind));
    
    %%%%%%%%%%% Flag HL %%%%%%%%%%%
    % Hot load check flag
    if (abs(calibratedSpectra(i).THot - calibrationTool.THotTh)>calibrationTool.THotAbsThresh | calibratedSpectra(i).stdTHot>calibrationTool.hotTemperatureStdThreshold)
        hotLoadOK=0;
    else
        hotLoadOK=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cold and Hot spectra variation among this cycle
    calibratedSpectra(i).meanStdHotSpectra=nanmean(calibratedSpectra(i).stdHotSpectra);
    calibratedSpectra(i).meanStdColdSpectra=nanmean(calibratedSpectra(i).stdColdSpectra);
    
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
    
    %%%%%%%%%%% Flag Pointing %%%%%%%%%%%
    % Antenna angle
    calibratedSpectra(i).meanAngleAntenna=mean(logFile.Elevation_Angle(ia));
    calibratedSpectra(i).stdAngleAntenna=std(logFile.Elevation_Angle(ia));
    

    if calibratedSpectra(i).stdAngleAntenna > calibrationTool.stdAntAngleThresh
        PointingAngleOK=0;
    else
        PointingAngleOK=1;
    end
    
%     %%%%%%%%%%% Flag 9 %%%%%%%%%%%
%     % FFTS aquisition variable
%     calibratedSpectra(i).FFT_adc_range=logFile.FFT_adc_range(1);
%     if sum(logFile.FFT_adc_overload(ind))>0
%         FFT_adc_overload_OK=0;
%     else
%         FFT_adc_overload_OK=1;
%     end
    
    % Other variables (if existing)
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
    % Correspond to the first sky measurements taken into account for the
    % mean calibrated spectra.
    
    % As we are always using daily raw files:
    calibratedSpectra(i).year=calibratedSpectra(i).theoreticalStartTime.Year;
    calibratedSpectra(i).month=calibratedSpectra(i).theoreticalStartTime.Month;
    calibratedSpectra(i).day=calibratedSpectra(i).theoreticalStartTime.Day;
    
    calibratedSpectra(i).timeMin = datestr(calibratedSpectra(i).theoreticalStartTime,'YYYY_mm_dd_HH:MM:SS');
    calibratedSpectra(i).timeMax = datestr(calibratedSpectra(i).theoreticalStartTime + minutes(calibratedSpectra(i).calibrationTime),'YYYY_mm_dd_HH:MM:SS');
    
    calibratedSpectra(i).timeMin=datenum(calibratedSpectra(i).timeMin,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    calibratedSpectra(i).timeMax=datenum(calibratedSpectra(i).timeMax,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    
    if ~isempty(ia)
        
        %only possible if ia is not empty !
        calibratedSpectra(i).firstSkyTime=logFile.time(ia(1))-datenum(1970,1,1);
        calibratedSpectra(i).lastSkyTime=logFile.time(ia(end))-datenum(1970,1,1);
        
        % "mean time" of the calibration cycle (mean of all antenna measurements)
        calibratedSpectra(i).meanAntTime = nanmean(logFile.dateTime(ia));
        calibratedSpectra(i).meanDatetime = datenum(calibratedSpectra(i).meanAntTime)-datenum(1970,1,1);
        
        calibratedSpectra(i).timeOfDay = 24*(datenum(calibratedSpectra(i).meanAntTime) -datenum(calibratedSpectra(i).year,calibratedSpectra(i).month,calibratedSpectra(i).day));

%     if calibratedSpectra(i).theoreticalStartTime.Month < 10
%         m = ['0' num2str(calibratedSpectra(i).theoreticalStartTime.Month)];
%     else
%         m = num2str(calibratedSpectra(i).theoreticalStartTime.Month);
%     end
%     
%     if calibratedSpectra(i).theoreticalStartTime.Day < 10
%         d = ['0' num2str(calibratedSpectra(i).theoreticalStartTime.Day)];
%     else
%         d = num2str(calibratedSpectra(i).theoreticalStartTime.Day);
%     end
%     
%     calibratedSpectra(i).date=[num2str(logFile.Year(1)) '_' m '_' d];
    
    %calibratedSpectra(i).startTimeInt8=int8(datestr(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ'));
    
    else
        calibratedSpectra(i).firstSkyTime = -9999;
        calibratedSpectra(i).lastSkyTime = -9999;
        calibratedSpectra(i).timeOfDay = -9999;
        
        calibratedSpectra(i).meanAntTime = calibratedSpectra(i).theoreticalStartTime + 0.5*minutes(calibratedSpectra(i).calibrationTime);
        calibratedSpectra(i).meanDatetime = datenum(calibratedSpectra(i).meanAntTime)-datenum(1970,1,1);
    end
      
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
        systemTemperatureOK,...
        LN2SensorsOK,...
        LN2LevelOK,...
        hotLoadOK,...
        PointingAngleOK];
    
    % Error vector description:
    calibratedSpectra(i).errorVectorDescription=[
        "sufficientNumberOfIndices",...
        "systemTemperatureOK",...
        "LN2SensorsOK",...
        "LN2LevelOK",...
        "hotLoadOK",...
        "PointingAngleOK"];
        
    calibratedSpectra(i).outlierCalib = NaN;
    warning('off','backtrace')
    if (sum(calibratedSpectra(i).errorVector)<length(calibratedSpectra(i).errorVector))
        calibratedSpectra(i).outlierCalib = 1;
        errorV=num2str(calibratedSpectra(i).errorVector);
        disp(['Calibration Cycle number ' num2str(i) ', TOD: ' datestr(timeofday(calibratedSpectra(i).meanAntTime),'HH:MM:SS')])
        warning(['Problem with this calibration, error code : ' errorV]);
        disp(calibratedSpectra(i).errorVectorDescription(~calibratedSpectra(i).errorVector))
    end
end


%calibratedSpectra(emptyTimestamp)=[];

% Ording by datenum
%T=struct2table(calibratedSpectra);
%sortedT=sortrows(T,'meanDatetime');
%calibratedSpectra=table2struct(sortedT);

% figure
% clf
% set(gcf, 'PaperPosition', [1 1 19 27.7])
% subplot(3,2,1);
% plot([calibratedSpectra.globalHotCounts],'r');
% hold on
% plot([calibratedSpectra.globalColdCounts],'b');
% plot([calibratedSpectra.globalAntennaCounts],'g');
% ylabel('Counts []') 
% subplot(3,2,2);
% plot([calibratedSpectra.globalTa],'g');
% ylabel('Ta [K]');
% subplot(3,2,3);
% plot([calibratedSpectra.globalTN],'k');
% ylabel('TN [K]');
% subplot(3,2,4);
% plot([calibratedSpectra.allTHot],'r');
% ylabel('THot [K]');
% ylim([290,300])
% subplot(3,2,5);
% for i = 1:length(calibratedSpectra)
%     %plot([calibratedSpectra(i).freq],calibratedSpectra(i).Tb);
%     plot(calibratedSpectra(i).Tb);
%     hold on
%     %xlim([0,1000])
%     ylim([50,250])
%     ylabel('T_B [K]');
% end   
% subplot(3,2,6);
% for i = 1:length(calibratedSpectra)
%     %plot([calibratedSpectra(i).freq],calibratedSpectra(i).TSys);
%     plot(calibratedSpectra(i).TSys);
%     hold on
%     ylabel('TN [K]');
%     ylim([450,600])
%     %xlim([0,1000])
% end
% print([calibrationTool.level1Folder calibrationTool.dateStr '_' calibrationTool.spectrometer],'-dpdf','-fillpage')
% close
end

