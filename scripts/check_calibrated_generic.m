function calibratedSpectra = check_calibrated_generic(log,retrievalTool,calibratedSpectra)
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
% General warning raised during the calibration
for i = 1:size(calibratedSpectra,2)
    ia=calibratedSpectra(i).antennaInd;
    ih=calibratedSpectra(i).hotInd;
    ic=calibratedSpectra(i).coldInd;
    
    ind=sort([ia ih ic]);
    
    % The number of indices for the 3 positions considered in the
    % calibration cycle should be more than a certain threshold
    if ((length(ih)>retrievalTool.minNumberOfIndicePerCycle) || (length(ia) > retrievalTool.minNumberOfIndicePerCycle) || (length(ic)>retrievalTool.minNumberOfIndicePerCycle))
        calibratedSpectra(i).sufficientNumberOfIndices=1;
    else
        calibratedSpectra(i).sufficientNumberOfIndices=0;
        warning('number of spectra low');
    end
    
    % Effective calibration time for this cycle (TO CHECK IF NEEDED ?)
    effectiveTime=0;
    %if i<size(calibratedSpectra,2)
    for l= 1:length(ind)-1
        effectiveTime=effectiveTime+(log.t(ind(l+1))-log.t(ind(l)))*60;
    end
    % adding the last time if possible
    if i<size(calibratedSpectra,2)
        effectiveTime=effectiveTime+(log.t(ind(l+1)+1)-log.t(ind(l+1)));
    else
        % TODO
        effectiveTime=effectiveTime;
    end
    calibratedSpectra(i).effectiveCalibrationTime=effectiveTime;
      
    %calibratedSpectra(i).effectiveCalibrationTimeHot=length(ih)*retrievalTool.calibTimeHot;
    %calibratedSpectra(i).effectiveCalibrationTimeAntenna=length(ia)*retrievalTool.calibTimeAntenna;
    %calibratedSpectra(i).effectiveCalibrationTimeCold=length(ic)*retrievalTool.calibTimeCold;
    %calibratedSpectra(i).effectiveCalibrationTime=calibratedSpectra(i).effectiveCalibrationTimeHot+calibratedSpectra(i).effectiveCalibrationTimeAntenna+calibratedSpectra(i).effectiveCalibrationTimeCold;
    
    % Frequency vector
    calibratedSpectra(i).if = retrievalTool.samplingRateFFTS/2 * [0:1/retrievalTool.numberOfChannels:1-1/retrievalTool.numberOfChannels];
    
    calibratedSpectra(i).observationFreq=retrievalTool.observationFreq;
        
    calibratedSpectra(i).LOFreqTot=retrievalTool.LOFreqTot;
    
    calibratedSpectra(i).freq=calibratedSpectra(i).if*1e6+calibratedSpectra(i).LOFreqTot;
    
    calibratedSpectra(i).df=retrievalTool.samplingRateFFTS/(2*retrievalTool.numberOfChannels);
    
    
    %bw=retrievalTool.instrumentBandwidth;
    %nChannel=retrievalTool.numberOfChannels;
    %df=bw/(nChannel+1); % TOCHECK
    % lc=log.Spectr_line_center(1);
    %calibratedSpectra(i).freq=horzcat(sort(calibratedSpectra(i).LOFreqTot-df*(0:retrievalTool.DCChannel-1)),calibratedSpectra(i).LOFreqTot+df*(1:nChannel-retrievalTool.DCChannel));
    %calibratedSpectra(i).if=calibratedSpectra(i).freq-calibratedSpectra(i).freq(1);
    %calibratedSpectra(i).freq=(calibratedSpectra(i).f0-(lc*df)):df:calibratedSpectra(i).f0+((nChannel-(lc+1))*df);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Antenna angle
    calibratedSpectra(i).meanAngleAntenna=mean(log.Elevation_Angle(ia));
    calibratedSpectra(i).stdAngleAntenna=std(log.Elevation_Angle(ia));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hot load check flag
    if (abs(calibratedSpectra(i).THot - retrievalTool.THotTh)>retrievalTool.THotAbsThresh | calibratedSpectra(i).stdTHot>retrievalTool.hotTemperatureStdThreshold)
        calibratedSpectra(i).hotLoadOK=0;
    else
        calibratedSpectra(i).hotLoadOK=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % System Temperature
    % Computing TN around the line center (approximately +- 200 MHz)
    centerChannels=find(calibratedSpectra(i).freq>=calibratedSpectra(i).observationFreq-200e6 & calibratedSpectra(i).freq<calibratedSpectra(i).observationFreq+200e6);
    
    Ycenter=calibratedSpectra(i).Yspectral(centerChannels);
    
    % Removing extreme outliers before computing TN:
    boxCarFilter=ones(100,1)/100;
    Yfiltered=conv(Ycenter,boxCarFilter,'same');
    YcenterSmoothed=Ycenter(abs(Ycenter-Yfiltered)<2*nanstd(Ycenter));
    
    TSysCenter=(calibratedSpectra(i).THot-YcenterSmoothed*retrievalTool.TCold)./(YcenterSmoothed-1);
    
    %TSysCenter2=(calibratedSpectra(i).THot-mean(Ycenter)*retrievalTool.TCold)./(mean(Ycenter)-1);
    
    calibratedSpectra(i).TSys=nanmean(TSysCenter);

    if (abs(calibratedSpectra(i).TSys-retrievalTool.TSysCenterTh)>retrievalTool.TSysThresh | calibratedSpectra(i).stdTSys > retrievalTool.stdTSysThresh)
        calibratedSpectra(i).systemTemperatureOK=0;
    else
        calibratedSpectra(i).systemTemperatureOK=1;
    end
    
    % Tsys from the log
    calibratedSpectra(i).TSysLog=nanmean(log.FE_T_Sys(ind));
    calibratedSpectra(i).stdTSysLog=nanstd(log.FE_T_Sys(ind));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cold and Hot spectra variation among this cycle
    calibratedSpectra(i).meanStdHotSpectra=nanmean(calibratedSpectra(i).stdHotSpectra);
    calibratedSpectra(i).meanStdColdSpectra=nanmean(calibratedSpectra(i).stdColdSpectra);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Liquid Nitrogen
    if ~all(log.LN2_Sensors_OK(ind)==1)
        calibratedSpectra(i).LN2SensorsOK=0;
    else
        calibratedSpectra(i).LN2SensorsOK=1;
    end
    if ~all(log.LN2_Level_OK(ind))
        calibratedSpectra(i).LN2LevelOK=0;
    else
        calibratedSpectra(i).LN2LevelOK=1;
    end
    
    % Other variables (if existing)
    if isfield(log,'T_Room')
        calibratedSpectra(i).TempRoom=nanmean(log.T_Room(ind));
        calibratedSpectra(i).stdTempRoom=nanstd(log.T_Room(ind));
    else
        calibratedSpectra(i).TempRoom=-9999;
        calibratedSpectra(i).stdTempRoom=-9999;
    end
    if isfield(log,'T_Out')
        calibratedSpectra(i).TempOut=nanmean(log.T_Out(ind));
        calibratedSpectra(i).stdTempOut=nanstd(log.T_Out(ind));
    else
        calibratedSpectra(i).TempOut=-9999;
        calibratedSpectra(i).stdTempOut=-9999;
    end
    if isfield(log,'T_Window')
        calibratedSpectra(i).TempWindow=nanmean(log.T_Window(ind));
        calibratedSpectra(i).stdTempWindow=nanstd(log.T_Window(ind));
    else
        calibratedSpectra(i).TempWindow=-9999;
        calibratedSpectra(i).stdTempWindow=-9999;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FFTS aquisition variable
    calibratedSpectra(i).FFT_adc_range=log.FFT_adc_range(1);
    if sum(log.FFT_adc_overload(ind))>0
        calibratedSpectra(i).FFT_adc_overload_OK=0;
    else
        calibratedSpectra(i).FFT_adc_overload_OK=1;
    end
    
    % FFTS number of aquisition
    if ((all(log.FFT_Nr_of_acq(ia))==retrievalTool.numberOfAquisitionSpectraAntenna) & (all(log.FFT_Nr_of_acq(ih))==retrievalTool.numberOfAquisitionSpectraHot) & (all(log.FFT_Nr_of_acq(ih))==retrievalTool.numberOfAquisitionSpectraCold))
        calibratedSpectra(i).FFT_nr_aquisition_OK=0;
    else
        calibratedSpectra(i).FFT_nr_aquisition_OK=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the start and stop time for this calibration cycle
    % Correspond to the first sky measurements taken into account for the
    % mean calibrated spectra.
       
    calibratedSpectra(i).dateStart=datestr(log.x(1:6,ia(1))','yyyymmddTHHMMSSZ');
    calibratedSpectra(i).dateStop=datestr(log.x(1:6,ia(end))','yyyymmddTHHMMSSZ');
    
    calibratedSpectra(i).firstSkyTime=datenum(calibratedSpectra(i).dateStart,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
    calibratedSpectra(i).lastSkyTime=datenum(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ')-datenum(1970,1,1);
    
    % As we are always using daily raw files:
    calibratedSpectra(i).year=log.Year(1);
    calibratedSpectra(i).month=log.Month(1);
    calibratedSpectra(i).day=log.Day(1);
    
    if log.Month(1) < 10
        m = ['0' num2str(log.Month(1))];
    else
        m = num2str(log.Month(1));
    end
    
    if log.Day(1) < 10
        d = ['0' num2str(log.Day(1))];
    else
        d = num2str(log.Day(1));
    end
    
    calibratedSpectra(i).date=[num2str(log.Year(1)) '_' m '_' d];
    
    % as well as the "mean time" of the calibration cycle (mean of all
    % antenna measurements)
    meanDatetime=[calibratedSpectra(i).date '_' datestr(mean(log.t(ia))/24,'HH:MM:SS')];
    
    theoreticalminTime=[calibratedSpectra(i).date '_' datestr(calibratedSpectra(i).theoreticalStartTime/24,'HH:MM:SS')];
    %theoreticalmaxTime=[calibratedSpectra(i).date '_' datestr((calibratedSpectra(i).theoreticalStartTime+calibratedSpectra(i).calibrationTime/60)/24,'HH:MM:SS')];
    
    calibratedSpectra(i).timeMin=datenum(theoreticalminTime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    
    calibratedSpectra(i).meanDatetime=datenum(meanDatetime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
    calibratedSpectra(i).meanDatetimeUnit='days since 1970-01-01 00:00:00';
    calibratedSpectra(i).calendar='standard';
    calibratedSpectra(i).timeOfDay=mean(log.t(ia));
    
    %calibratedSpectra(i).startTimeInt8=int8(datestr(calibratedSpectra(i).dateStop,'yyyymmddTHHMMSSZ'));
   
    
    % Error vector for this calibration cycle
    calibratedSpectra(i).errorVector=[calibratedSpectra(i).sufficientNumberOfIndices,...
        calibratedSpectra(i).systemTemperatureOK,...
        calibratedSpectra(i).hotAngleRemoved,...
        calibratedSpectra(i).coldAngleRemoved,...
        calibratedSpectra(i).antennaAngleRemoved,...
        calibratedSpectra(i).LN2SensorsOK,...
        calibratedSpectra(i).LN2LevelOK,...
        calibratedSpectra(i).hotLoadOK,...
        calibratedSpectra(i).FFT_adc_overload_OK];
    
    if (sum(calibratedSpectra(i).errorVector([1,2,6,7,8,9]))<6 | sum(calibratedSpectra(i).errorVector([3,4,5]))>0)
        errorV=num2str(calibratedSpectra(i).errorVector);
        disp(['Calibration Cycle number ' num2str(i) ', TOD: ' num2str(calibratedSpectra(i).timeOfDay)])
        warning(['Problem with this calibration, error code : ' errorV]);
    end
end

end

