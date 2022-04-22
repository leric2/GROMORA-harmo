function [o, calibratedSpectra] = run_balancing_calibration(rawSpectra,logFile,calibrationTool,~)
o = [];

disp(['run balancing calibration for ' calibrationTool.dateStr])

if length(unique(logFile.dir))>1
    1
    calType = 'left';
    disp(['Calibrate ' calType])
    S1 = balancing_calibration_generic(rawSpectra,logFile,calibrationTool,calType);
    2
    calType = 'right';
    disp(['Calibrate ' calType])
    S2 = balancing_calibration_generic(rawSpectra,logFile,calibrationTool,calType);
    
    calibratedSpectra = [S1 S2];
else
    0
    calType = '0';
    calibratedSpectra = balancing_calibration_generic(rawSpectra,logFile,calibrationTool,calType);
end
    
















