function calibratedSpectra = check_calibrated_generic(standardLog,calibrationTool,calibratedSpectra)
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
% ARGUMENTS | INPUTS:   - standardLog: harmonized GROSOM log file 
%           |           - calibrationTool
%           |           - calibratedSpectra
%           |
%           |
%           | OUTPUTS: - calibratedSpectra
%           |
%           |
% CALLS     | 
%           | 
%           | 
%           | 
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
    if ((length(ih)>calibrationTool.minNumberOfIndicePerCycle) || (length(ia) > calibrationTool.minNumberOfIndicePerCycle) || (length(ic)>calibrationTool.minNumberOfIndicePerCycle))
        sufficientNumberOfIndices=1;
    else
        sufficientNumberOfIndices=0;
        warning('Low number of spectra for this cycle');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency vector
    calibratedSpectra(i).if = calibrationTool.samplingRateFFTS/2 * [0:1/calibrationTool.numberOfChannels:1-1/calibrationTool.numberOfChannels];
    
    calibratedSpectra(i).observationFreq=calibrationTool.observationFreq;
        
    calibratedSpectra(i).LOFreqTot=calibrationTool.LOFreqTot;
    
    calibratedSpectra(i).freq=calibratedSpectra(i).if*1e6+calibratedSpectra(i).LOFreqTot;
    
    calibratedSpectra(i).df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    
    %bw=retrievalTool.instrumentBandwidth;
    %nChannel=retrievalTool.numberOfChannels;
    %df=bw/(nChannel+1); % TOCHECK
    % lc=log.Spectr_line_center(1);
    %calibratedSpectra(i).freq=horzcat(sort(calibratedSpectra(i).LOFreqTot-df*(0:retrievalTool.DCChannel-1)),calibratedSpectra(i).LOFreqTot+df*(1:nChannel-retrievalTool.DCChannel));
    %calibratedSpectra(i).if=calibratedSpectra(i).freq-calibratedSpectra(i).freq(1);
    %calibratedSpectra(i).freq=(calibratedSpectra(i).f0-(lc*df)):df:calibratedSpectra(i).f0+((nChannel-(lc+1))*df);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Antenna angle
    calibratedSpectra(i).meanAngleAntenna=mean(standardLog.Elevation_Angle(ia));
    calibratedSpectra(i).stdAngleAntenna=std(standardLog.Elevation_Angle(ia));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % System Temperature
    % Computing TN around the line center (approximately +- 200 MHz)
    centerChannels=find(calibratedSpectra(i).freq>=calibratedSpectra(i).observationFreq-200e6 & calibratedSpectra(i).freq<calibratedSpectra(i).observationFreq+200e6);
    
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
    calibratedSpectra(i).TSysLog=nanmean(standardLog.FE_T_Sys(ind));
    calibratedSpectra(i).stdTSysLog=nanstd(standardLog.FE_T_Sys(ind));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cold and Hot spectra variation among this cycle
    calibratedSpectra(i).meanStdHotSpectra=nanmean(calibratedSpectra(i).stdHotSpectra);
    calibratedSpectra(i).meanStdColdSpectra=nanmean(calibratedSpectra(i).stdColdSpectra);
    
    %%%%%%%%%%% Flag 6 %%%%%%%%%%%    
    % Liquid Nitrogen sensors
    if ~all(standardLog.LN2_Sensors_OK(ind)==1)
        LN2SensorsOK=0;
    else
        LN2SensorsOK=1;
    end
    
    %%%%%%%%%%% Flag 6 %%%%%%%%%%%
    % Liquid Nitrogen level
    if ~all(standardLog.LN2_Level_OK(ind))
        LN2LevelOK=0;
    else
        LN2LevelOK=1;
    end

    %%%%%%%%%%% Flag 8 %%%%%%%%%%%
    % Hot load check flag
    if (abs(calibratedSpectra(i).THot - calibrationTool.THotTh)>calibrationTool.THotAbsThresh | calibratedSpectra(i).stdTHot>calibrationTool.hotTemperatureStdThreshold)
        hotLoadOK=0;
    else
        hotLoadOK=1;
    end
    
    %%%%%%%%%%% Flag 9 %%%%%%%%%%%
    % FFTS aquisition variable
    calibratedSpectra(i).FFT_adc_range=standardLog.FFT_adc_range(1);
    if sum(standardLog.FFT_adc_overload(ind))>0
        FFT_adc_overload_OK=0;
    else
        FFT_adc_overload_OK=1;
    end
    
    %%%%%%%%%%% Flag ... %%%%%%%%%%%    NOT USED 
    % FFTS number of aquisition
    if ((all(standardLog.FFT_Nr_of_acq(ia))==calibrationTool.numberOfAquisitionSpectraAntenna) & (all(standardLog.FFT_Nr_of_acq(ih))==calibrationTool.numberOfAquisitionSpectraHot) & (all(standardLog.FFT_Nr_of_acq(ih))==calibrationTool.numberOfAquisitionSpectraCold))
        calibratedSpectra(i).FFT_nr_aquisition_OK=1;
    else
        calibratedSpectra(i).FFT_nr_aquisition_OK=0;
    end
    
    % Error vector for this calibration cycle
    calibratedSpectra(i).errorVector=[
        sufficientNumberOfIndices,...
        systemTemperatureOK,...
        length(ih),...
        length(ic),...
        length(ia),...
        LN2SensorsOK,...
        LN2LevelOK,...
        hotLoadOK,...
        FFT_adc_overload_OK];
    
    % Error vector description:
    calibratedSpectra(i).errorVectorDescription=[
        "sufficientNumberOfIndices",...
        "systemTemperatureOK",...
        "numberHotIndice",...
        "numberColdIndice",...
        "numberAntennaIndice",...
        "LN2SensorsOK",...
        "LN2LevelOK",...
        "hotLoadOK",...
        "FFT_adc_overload_OK"];
    
    if (sum(calibratedSpectra(i).errorVector([1,2,6,7,8,9]))<6 | sum(calibratedSpectra(i).errorVector([3,4,5]))>0)
        errorV=num2str(calibratedSpectra(i).errorVector);
        disp(['Calibration Cycle number ' num2str(i) ', TOD: ' num2str(calibratedSpectra(i).timeOfDay)])
        warning(['Problem with this calibration, error code : ' errorV]);
        disp(calibratedSpectra(i).errorVectorDescription)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding other variables (if existing)
    if isfield(standardLog,'T_Room')
        calibratedSpectra(i).TempRoom=nanmean(standardLog.T_Room(ind));
        calibratedSpectra(i).stdTempRoom=nanstd(standardLog.T_Room(ind));
    else
        calibratedSpectra(i).TempRoom=-9999;
        calibratedSpectra(i).stdTempRoom=-9999;
    end
    if isfield(standardLog,'T_Out')
        calibratedSpectra(i).TempOut=nanmean(standardLog.T_Out(ind));
        calibratedSpectra(i).stdTempOut=nanstd(standardLog.T_Out(ind));
    else
        calibratedSpectra(i).TempOut=-9999;
        calibratedSpectra(i).stdTempOut=-9999;
    end
    if isfield(standardLog,'T_Window')
        calibratedSpectra(i).TempWindow=nanmean(standardLog.T_Window(ind));
        calibratedSpectra(i).stdTempWindow=nanstd(standardLog.T_Window(ind));
    else
        calibratedSpectra(i).TempWindow=-9999;
        calibratedSpectra(i).stdTempWindow=-9999;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time variable for this cycle
    % Correspond to the first sky measurements taken into account for the
    % mean calibrated spectra.
       
    calibratedSpectra(i).dateStart=datestr(standardLog.x(1:6,ia(1))','yyyymmddTHHMMSSZ');
    calibratedSpectra(i).dateStop=datestr(standardLog.x(1:6,ia(end))','yyyymmddTHHMMSSZ');
    
    calibratedSpectra(i).firstSkyTime=datenum(calibratedSpectra(i).dateStart,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
    calibratedSpectra(i).lastSkyTime=datenum(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
    
    % As we are always using daily raw files:
    calibratedSpectra(i).year=standardLog.Year(1);
    calibratedSpectra(i).month=standardLog.Month(1);
    calibratedSpectra(i).day=standardLog.Day(1);
    
    if standardLog.Month(1) < 10
        m = ['0' num2str(standardLog.Month(1))];
    else
        m = num2str(standardLog.Month(1));
    end
    
    if standardLog.Day(1) < 10
        d = ['0' num2str(standardLog.Day(1))];
    else
        d = num2str(standardLog.Day(1));
    end
    
    calibratedSpectra(i).date=[num2str(standardLog.Year(1)) '_' m '_' d];
    
    % as well as the "mean time" of the calibration cycle (mean of all
    % antenna measurements)
    meanDatetime=[calibratedSpectra(i).date '_' datestr(mean(standardLog.t(ia))/24,'HH:MM:SS')];
    
    theoreticalminTime=[calibratedSpectra(i).date '_' datestr(calibratedSpectra(i).theoreticalStartTime/24,'HH:MM:SS')];
    %theoreticalmaxTime=[calibratedSpectra(i).date '_' datestr((calibratedSpectra(i).theoreticalStartTime+calibratedSpectra(i).calibrationTime/60)/24,'HH:MM:SS')];
    
    calibratedSpectra(i).timeMin=datenum(theoreticalminTime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    
    calibratedSpectra(i).meanDatetime=datenum(meanDatetime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    %calibratedSpectra(i).meanDatetimeUnit='days since 1970-01-01 00:00:00';
    %calibratedSpectra(i).calendar='standard';
    calibratedSpectra(i).timeOfDay=mean(standardLog.t(ia));
    
    %calibratedSpectra(i).startTimeInt8=int8(datestr(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ'));
   
    
        
    % Effective calibration time for this cycle (TO CHECK IF NEEDED ?)
    effectiveTime=0;
    %if i<size(calibratedSpectra,2)
    for l= 1:length(ind)-1
        effectiveTime=effectiveTime+(standardLog.t(ind(l+1))-standardLog.t(ind(l)))*60;
    end
    % adding the last time if possible
    if i<size(calibratedSpectra,2)
        effectiveTime=effectiveTime+(standardLog.t(ind(l+1)+1)-standardLog.t(ind(l+1)));
    else
        % TODO
        effectiveTime=effectiveTime;
    end
    calibratedSpectra(i).effectiveCalibrationTime=effectiveTime;
      
    %calibratedSpectra(i).effectiveCalibrationTimeHot=length(ih)*retrievalTool.calibTimeHot;
    %calibratedSpectra(i).effectiveCalibrationTimeAntenna=length(ia)*retrievalTool.calibTimeAntenna;
    %calibratedSpectra(i).effectiveCalibrationTimeCold=length(ic)*retrievalTool.calibTimeCold;
    %calibratedSpectra(i).effectiveCalibrationTime=calibratedSpectra(i).effectiveCalibrationTimeHot+calibratedSpectra(i).effectiveCalibrationTimeAntenna+calibratedSpectra(i).effectiveCalibrationTimeCold;
   
    
    
    
   
end

end

