function [drift,calibratedSpectra] = calibrate_balancing_with_checks(rawSpectra,logFile,calibrationTool,calType)
% NAME          | calibrate_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat adapted by Alistair Bell
% CREATION      | 01.2023
%               |
% ABSTRACT      | Doing a hot-cold calibration with additional checks based
%               | on the function calibration_generic.
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: 1. rawSpectra: Matrix of raw data (1 line per
%               |           cycle, channels as columns.
%               |         2. logFile: standardized log file
%               |         3. calibrationTool:
%               |               - calibrationVersion
%               |               - calibrationTime
%               |               - instrumentName
%               |               - indiceHot, indiceAntenna, indiceCold
%               |               - TCold
%               |               - Day, Month, Year
%               |               - timeZone
%               |               - elevationAngleHot, elevationAngleCold
%               |               - elevationAngleHotTol
%               |               - elevationAngleColdTol
%               |               - threshNumRawSpectraHot
%               |               - threshNumRawSpectraCold
%               |               - elevationAngleAntenna
%               |               - elevationAngleTolerance
%               |               - threshNumRawSpectraAnt
%               |               - numberOfChannels
%               |               - outlierDectectionType
%               |              
%               |         4. calType: Calibration type to do
%               |               1. standard: mean antenna vs mean hot/cold 
%               |                  spectra
%               |               2. debug: perform the standard one and:
%               |                  - Mean Up/Down: mean AT up, mean AT Down
%               |                  vs mean hot/cold
%               |                  - all cycles: all individual
%               |                  hot-cold-sky cycle BUT uncleaned for 
%               |                  potential spurious spectra ! Caution
%               |                  when comparing it to 'standard' mode if
%               |                  outliers in hot-cold are identified.
%               |                  
%               |
%               | OUTPUTS: 1. drift: structure containing mean values of all
%               |               channels for this day
%               |          2. calibratedSpectra: "structure array" of
%               |               calibrated data
%               |
% CALLS         | find_up_down_cycle()
%               |
%==========================================================================

% Calibration version
calibVersion = calibrationTool.calibrationVersion;

% CalibrationTime in Minute
calibTime=calibrationTool.calibrationTime;

%===== initialize output
disp('starting cold countstokelvin')
T_rec=nan(size(data.S_line));

%===== load physical constants
%  physical_constants
RE = 6371000;

% Extract all indices corresponding to hot, cold and sky observations in
% the log file.
initialIndices={
    find(logFile.Position==calibrationTool.indiceHot & logFile.Tipping_Curve_active==0);      % Hot
    find(logFile.Position==calibrationTool.indiceAntenna & logFile.Tipping_Curve_active==0);  % Antenna
    find(logFile.Position==calibrationTool.indiceCold & logFile.Tipping_Curve_active==0);     % Cold
    find(logFile.Position==calibrationTool.indiceCold & logFile.Tipping_Curve_active==0);     % Reference
    };
    
driftIndices = initialIndices;
standardOrder = cn(1)/cn(2)>0.9 && cn(3)/cn(2)>0.9;
if ~standardOrder
    warning('order of the cycle is not standard, no drift quantities !')
else
    % we remove extra indices to keep only the smallest number of indices
    %between hot, cold and sky
    for i=1:length(initialIndices); 	driftIndices{i}=initialIndices{i}(1:min(cn)); end; 
end



















