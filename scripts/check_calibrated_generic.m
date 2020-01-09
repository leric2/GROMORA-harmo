function calibratedSpectra = check_calibrated_generic(log,retrievalTool,calibratedSpectra)
%CHECK_CALIBRATED_SPECTRA Quality check of the calibrated spectra
%
% Completing and quality checking each calibration cycle

for i = 1:size(calibratedSpectra,2)
    ia=calibratedSpectra(i).antennaInd;
    ih=calibratedSpectra(i).hotInd;
    ic=calibratedSpectra(i).coldInd;
    
    ind=sort([ia; ih; ic]);
    
    % The number of indices for the 3 cases should be the same
    if ~(length(ih) == length(ia) == length(ic))
        calibratedSpectra(i).NumberOfIndicesForEachCaseFlag=1;
    else
        calibratedSpectra(i).NumberOfIndicesForEachCaseFlag=0;
    end
    
    % Antenna angle check
    calibratedSpectra(i).meanAngleAntenna=mean(log.Elevation_Angle(ia));
    calibratedSpectra(i).stdAngleAntenna=std(log.Elevation_Angle(ia));
    
    if ((calibratedSpectra(i).meanAngleAntenna>retrievalTool.elevationAngleAntenna-retrievalTool.elevationAngleTolerance) || (calibratedSpectra(i).meanAngleAntenna<retrievalTool.elevationAngleAntenna+retrievalTool.elevationAngleTolerance))
        calibratedSpectra(i).angleFlag=1;
    else
        calibratedSpectra(i).angleFlag=0;
    end
    
    % Hot load check (already computed during calibration)
    if (calibratedSpectra(i).stdTHot>10)
        calibratedSpectra(i).hotLoadFlag=1;
    else
        calibratedSpectra(i).hotLoadFlag=0;
    end
        
    % Tsys
    calibratedSpectra(i).Tsys=nanmean(log.FE_T_Sys(ind));
    calibratedSpectra(i).stdTSys=nanstd(log.FE_T_Sys(ind));
    
    if (calibratedSpectra(i).stdTSys>10)
        calibratedSpectra(i).systemTemperatureFlag=1;
    else
        calibratedSpectra(i).systemTemperatureFlag=0;
    end
    
    % Liquid Nitrogen

    % Other variables ????
    calibratedSpectra(i).TempRoom=nanmean(log.AI_1(ind))*100;
    calibratedSpectra(i).TempOut=nanmean(log.AI_7(ind))*100;
    calibratedSpectra(i).TempWindow=nanmean(log.AI_3(ind))*100;
    
    calibratedSpectra(i).stdTempRoom=std(log.AI_1(ind)*100);
    calibratedSpectra(i).stdTempOut=std(log.AI_7(ind)*100);
    calibratedSpectra(i).stdTempWindow=std(log.AI_3(ind)*100);
    
    % saved the start and stop time for this calibration cycle
    calibratedSpectra(i).dateStart=log.x(1:6,ih(1))';
    calibratedSpectra(i).dateStop=log.x(1:6,ih(end))';
    
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
    
    % as well as the "mean time" of the day  /OR USE MEDIAN ?
    meanDatetime=[calibratedSpectra(i).date '_' datestr(mean(log.t(ih(1):ih(end)))/24,'HH:MM:SS')];
    calibratedSpectra(i).meanDatetime=datenum(meanDatetime,'YYYY_mm_dd_HH:MM:SS')-datenum(2000,1,1);

    calibratedSpectra(i).timeOfDay=mean(log.t(ih(1):ih(end)));
end

end

