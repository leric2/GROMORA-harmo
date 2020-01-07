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
    
    % Other variables ????
    calibratedSpectra(i).TempRoom=nanmean(log.AI_1(ind))*100;
    calibratedSpectra(i).TempOut=nanmean(log.AI_7(ind))*100;
    calibratedSpectra(i).TempWindow=nanmean(log.AI_3(ind))*100;
    
    calibratedSpectra(i).stdTempRoom=std(log.AI_1(ind)*100);
    calibratedSpectra(i).stdTempOut=std(log.AI_7(ind)*100);
    calibratedSpectra(i).stdTempWindow=std(log.AI_3(ind)*100);
    
    % saved the start and stop time for this calibration cycle
    calibratedSpectra(i).dateStart=log.x(1:6,ih(1));
    calibratedSpectra(i).dateStop=log.x(1:6,ih(end));
end

end

