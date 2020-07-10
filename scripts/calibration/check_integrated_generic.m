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
    ia=integratedSpectra(i).antennaIndCleanAngle;
    ih=integratedSpectra(i).hotInd;
    ic=integratedSpectra(i).coldInd;
    
    ind=sort([ia ih ic]);
    
    %%%%%%%%%%% Flag 1 %%%%%%%%%%%
    % The number of indices for the 3 positions:
    if ((length(ih) > calibrationTool.minNumberOfIndicePerCycle) &&...
            (length(ia) > calibrationTool.minNumberOfIndicePerCycle) &&...
            (length(ic) > calibrationTool.minNumberOfIndicePerCycle))
        sufficientNumberOfIndices=1;
    else
        sufficientNumberOfIndices=0;
        warning('Low number of spectra for this cycle');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency vector
    integratedSpectra(i).if = calibrationTool.samplingRateFFTS/2 * [0:1/calibrationTool.numberOfChannels:1-1/calibrationTool.numberOfChannels];
    
    integratedSpectra(i).observationFreq=calibrationTool.observationFreq;
        
    integratedSpectra(i).LOFreqTot=calibrationTool.LOFreqTot;
    
    integratedSpectra(i).freq=integratedSpectra(i).if*1e6+integratedSpectra(i).LOFreqTot;
    
    integratedSpectra(i).df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    
    %bw=retrievalTool.instrumentBandwidth;
    %nChannel=retrievalTool.numberOfChannels;
    %df=bw/(nChannel+1); % TOCHECK
    % lc=log.Spectr_line_center(1);
    %integratedSpectra(i).freq=horzcat(sort(integratedSpectra(i).LOFreqTot-df*(0:retrievalTool.DCChannel-1)),integratedSpectra(i).LOFreqTot+df*(1:nChannel-retrievalTool.DCChannel));
    %integratedSpectra(i).if=integratedSpectra(i).freq-integratedSpectra(i).freq(1);
    %integratedSpectra(i).freq=(integratedSpectra(i).f0-(lc*df)):df:integratedSpectra(i).f0+((nChannel-(lc+1))*df);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Antenna angle
    integratedSpectra(i).meanAngleAntenna=mean(standardLog.Elevation_Angle(ia));
    integratedSpectra(i).stdAngleAntenna=std(standardLog.Elevation_Angle(ia));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % System Temperature
    % Computing TN around the line center (approximately +- 200 MHz)
    centerChannels=find(integratedSpectra(i).freq>=integratedSpectra(i).observationFreq-200e6 & integratedSpectra(i).freq<integratedSpectra(i).observationFreq+200e6);
    
    Ycenter=integratedSpectra(i).Yspectral(centerChannels);
    
    % Removing extreme outliers before computing TN:
    boxCarFilter=ones(100,1)/100;
    Yfiltered=conv(Ycenter,boxCarFilter,'same');
    
    YcenterSmoothed=Ycenter(abs(Ycenter-Yfiltered)<2*nanstd(Ycenter));
    
    TSysCenter=(integratedSpectra(i).THot-YcenterSmoothed*calibrationTool.TCold)./(YcenterSmoothed-1);
    
    %TSysCenter2=(integratedSpectra(i).THot-mean(Ycenter)*retrievalTool.TCold)./(mean(Ycenter)-1);
    
    integratedSpectra(i).TSys=nanmean(TSysCenter);
    
    %%%%%%%%%%% Flag 2 %%%%%%%%%%%
    if (abs(integratedSpectra(i).TSys-calibrationTool.TSysCenterTh)>calibrationTool.TSysThresh | integratedSpectra(i).stdTSys > calibrationTool.stdTSysThresh)
        systemTemperatureOK=0;
    else
        systemTemperatureOK=1;
    end
    
    % Tsys from the log file
    integratedSpectra(i).TSysLog=nanmean(standardLog.FE_T_Sys(ind));
    integratedSpectra(i).stdTSysLog=nanstd(standardLog.FE_T_Sys(ind));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cold and Hot spectra variation among this cycle
    integratedSpectra(i).meanStdHotSpectra=nanmean(integratedSpectra(i).stdHotSpectra);
    integratedSpectra(i).meanStdColdSpectra=nanmean(integratedSpectra(i).stdColdSpectra);
    
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
    if (abs(integratedSpectra(i).THot - calibrationTool.THotTh)>calibrationTool.THotAbsThresh | integratedSpectra(i).stdTHot>calibrationTool.hotTemperatureStdThreshold)
        hotLoadOK=0;
    else
        hotLoadOK=1;
    end
    
    %%%%%%%%%%% Flag 9 %%%%%%%%%%%
    % FFTS aquisition variable
    integratedSpectra(i).FFT_adc_range=standardLog.FFT_adc_range(1);
    if sum(standardLog.FFT_adc_overload(ind))>0
        FFT_adc_overload_OK=0;
    else
        FFT_adc_overload_OK=1;
    end
    
    %%%%%%%%%%% Flag ... %%%%%%%%%%%    NOT USED 
    % FFTS number of aquisition
    if ((all(standardLog.FFT_Nr_of_acq(ia))==calibrationTool.numberOfAquisitionSpectraAntenna) & (all(standardLog.FFT_Nr_of_acq(ih))==calibrationTool.numberOfAquisitionSpectraHot) & (all(standardLog.FFT_Nr_of_acq(ih))==calibrationTool.numberOfAquisitionSpectraCold))
        integratedSpectra(i).FFT_nr_aquisition_OK=1;
    else
        integratedSpectra(i).FFT_nr_aquisition_OK=0;
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding other variables (if existing)
    if isfield(standardLog,'T_Room')
        integratedSpectra(i).TempRoom=nanmean(standardLog.T_Room(ind));
        integratedSpectra(i).stdTempRoom=nanstd(standardLog.T_Room(ind));
    else
        integratedSpectra(i).TempRoom=-9999;
        integratedSpectra(i).stdTempRoom=-9999;
    end
    if isfield(standardLog,'T_Out')
        integratedSpectra(i).TempOut=nanmean(standardLog.T_Out(ind));
        integratedSpectra(i).stdTempOut=nanstd(standardLog.T_Out(ind));
    else
        integratedSpectra(i).TempOut=-9999;
        integratedSpectra(i).stdTempOut=-9999;
    end
    if isfield(standardLog,'T_Window')
        integratedSpectra(i).TempWindow=nanmean(standardLog.T_Window(ind));
        integratedSpectra(i).stdTempWindow=nanstd(standardLog.T_Window(ind));
    else
        integratedSpectra(i).TempWindow=-9999;
        integratedSpectra(i).stdTempWindow=-9999;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time variable for this cycle
    % Correspond to the first sky measurements taken into account for the
    % mean calibrated spectra.
       
    integratedSpectra(i).dateStart=datestr(standardLog.x(1:6,ia(1))','yyyymmddTHHMMSSZ');
    integratedSpectra(i).dateStop=datestr(standardLog.x(1:6,ia(end))','yyyymmddTHHMMSSZ');
    
    integratedSpectra(i).firstSkyTime=datenum(integratedSpectra(i).dateStart,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
    integratedSpectra(i).lastSkyTime=datenum(integratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
    
    % As we are always using daily raw files:
    integratedSpectra(i).year=standardLog.Year(1);
    integratedSpectra(i).month=standardLog.Month(1);
    integratedSpectra(i).day=standardLog.Day(1);
    
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
    
    integratedSpectra(i).date=[num2str(standardLog.Year(1)) '_' m '_' d];
    
    % as well as the "mean time" of the calibration cycle (mean of all
    % antenna measurements)
    meanDatetime=[integratedSpectra(i).date '_' datestr(mean(standardLog.t(ia))/24,'HH:MM:SS')];
    
    theoreticalminTime=[integratedSpectra(i).date '_' datestr(integratedSpectra(i).theoreticalStartTime/24,'HH:MM:SS')];
    %theoreticalmaxTime=[integratedSpectra(i).date '_' datestr((integratedSpectra(i).theoreticalStartTime+integratedSpectra(i).calibrationTime/60)/24,'HH:MM:SS')];
    
    integratedSpectra(i).timeMin=datenum(theoreticalminTime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    
    integratedSpectra(i).meanDatetime=datenum(meanDatetime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    %integratedSpectra(i).meanDatetimeUnit='days since 1970-01-01 00:00:00';
    %integratedSpectra(i).calendar='standard';
    integratedSpectra(i).timeOfDay=mean(standardLog.t(ia));
    
    %integratedSpectra(i).startTimeInt8=int8(datestr(integratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ'));
    
    % Effective calibration time for this cycle (TO CHECK IF NEEDED ?)
    effectiveTime=0;
    %if i<size(integratedSpectra,2)
    for l= 1:length(ind)-1
        effectiveTime=effectiveTime+(standardLog.t(ind(l+1))-standardLog.t(ind(l)))*60;
    end
    % adding the last time if possible
    if i<size(integratedSpectra,2)
        effectiveTime=effectiveTime+(standardLog.t(ind(l+1)+1)-standardLog.t(ind(l+1)));
    else
        % TODO
        effectiveTime=effectiveTime;
    end
    integratedSpectra(i).effectiveCalibrationTime=effectiveTime;
      
    %integratedSpectra(i).effectiveCalibrationTimeHot=length(ih)*retrievalTool.calibTimeHot;
    %integratedSpectra(i).effectiveCalibrationTimeAntenna=length(ia)*retrievalTool.calibTimeAntenna;
    %integratedSpectra(i).effectiveCalibrationTimeCold=length(ic)*retrievalTool.calibTimeCold;
    %integratedSpectra(i).effectiveCalibrationTime=integratedSpectra(i).effectiveCalibrationTimeHot+integratedSpectra(i).effectiveCalibrationTimeAntenna+integratedSpectra(i).effectiveCalibrationTimeCold;
   
    integratedSpectra(i).numberOfIndices=[
        length(ih),...
        length(ic),...
        length(ia)];
    
    % Error vector for this calibration cycle
    integratedSpectra(i).errorVector=[
        sufficientNumberOfIndices,...
        systemTemperatureOK,...
        LN2SensorsOK,...
        LN2LevelOK,...
        hotLoadOK,...
        FFT_adc_overload_OK];
    
    % Error vector description:
    integratedSpectra(i).errorVectorDescription=[
        "sufficientNumberOfIndices",...
        "systemTemperatureOK",...
        "LN2SensorsOK",...
        "LN2LevelOK",...
        "hotLoadOK",...
        "FFT_adc_overload_OK"];
    
    if (sum(integratedSpectra(i).errorVector)<6)
        errorV=num2str(integratedSpectra(i).errorVector);
        disp(['Calibration Cycle number ' num2str(i) ', TOD: ' num2str(integratedSpectra(i).timeOfDay)])
        warning(['Problem with this calibration, error code : ' errorV]);
        disp(integratedSpectra(i).errorVectorDescription(~integratedSpectra(i).errorVector))
    end
   
end

end

