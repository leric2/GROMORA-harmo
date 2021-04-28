function calibrationTool = run_calibration(calibrationTool)
%==========================================================================
% NAME      | run_calibration.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 01.2020
%           |
% ABSTRACT  | The main function executing the calibration Level0 -> Level1a 
%           | for the instrument defined in calibrationTool. Some parts and
%           | functions are dependent on the instrument that we want to
%           | calibrate (also integrated in calibrationTool).
%           | 
%           |
% ARGUMENTS | INPUTS:  calibrationTool: structure containing all
%           | information about the calibration we want to perform.
%           | Documentation about this structure can be found in external
%           | document.
%           |
%           | OUTPUTS: calibrationTool with some additional fields.
%           |          
%           |
% SAVE      | level1a netCDF file and plots
%           |
% CALLS     | Some depends on instruments, all are stored in calibrationTool:
%           | - read_level0
%           | - harmonize_log
%           | - check_level0
%           | - reformat_spectra
%           | - flip_spectra
%           | - plot_raw_spectra
%           | - read_meteo_data
%           | - run_tipping_curve
%           | - calibrate
%           | - check_calibrated
%           | - plot_calibrated_spectra
%           | - save_level1a
%           |
% COMMENTS  | External documentation can be found for this function on the
%           | git server of IAP.
%==========================================================================
% Just checking that dateStr is a str...
assert(ischar(calibrationTool.dateStr),'Please enter the date in the right format')

% Check here that all required fields are filled in calibrationTool ?
% calibrationTool_complete(calibrationTool)

% Check if both bin and log file exist (does not check their content yet)
assert(exist([calibrationTool.file calibrationTool.logFileDataExtension],'file') ...
    && exist([calibrationTool.file calibrationTool.binaryDataExtension],'file'),'Raw files not found !')

% Start calibration
disp(['Starting the calibration process for ' calibrationTool.instrumentName ' ' calibrationTool.spectrometer ': ' calibrationTool.dateStr])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading level0 data...')

tic
[logFile,rawSpectra] = calibrationTool.read_level0(calibrationTool, 1);
toc

% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible (different function for each
% instrument).
logFile = calibrationTool.harmonize_log(calibrationTool, logFile);

% Reformat the raw spectra from vector to matrix
if size(rawSpectra,1) == 1
    rawSpectra = calibrationTool.reformat_spectra(rawSpectra, logFile, calibrationTool);
end

% Quality check of the raw data:
if calibrationTool.checkLevel0
    warningLevel0 = calibrationTool.check_level0(logFile, rawSpectra, calibrationTool);
else
    warningLevel0 = '';
end

% this function should be integrated to reformat_spectra()
if calibrationTool.flipped_spectra
    rawSpectra = calibrationTool.flip_spectra(rawSpectra, calibrationTool);
end

% Option for plotting spectra (to be improved...)
if calibrationTool.rawSpectraPlot
    calibrationTool.plot_raw_spectra(rawSpectra,0,1e4,20);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading meteo data and doing tipping curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logFile.meteo = calibrationTool.read_meteo_data(calibrationTool);

if calibrationTool.doTippingCurve
    % Tipping Curve
    TC = calibrationTool.run_tipping_curve(rawSpectra, logFile,...
        calibrationTool);
    if ~isempty(fieldnames(TC))
        logFile.TC = TC;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calibrating...')

% Calibration of the spectra using hot and cold load.  TODO: Add tipping
% There are different option for the calibration:
% - standard 
% - debug
% - time
% - all
[drift,calibratedSpectra] = calibrationTool.calibrate(rawSpectra, ...
    logFile, calibrationTool, calibrationTool.calType);

% Quality check of the calibrated spectra
% Also computing some additional metadata from the log file and storing
% everything in calibrated spectra
[calibratedSpectra, logFile] = calibrationTool.check_calibrated(logFile, ...
    calibrationTool, calibratedSpectra);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting and saving level 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calibrationTool.calibratedSpectraPlot
    try
        calibrationTool.plot_calibrated_spectra(calibrationTool, drift, ...
            logFile.meteo, calibratedSpectra, 24);
    catch ME
        warning(ME.identifier,'problem with the plotting:');
        disp(ME.message)
    end
end

% Saving calibrated spectra (level 1a) into NetCDF-4 file
if calibrationTool.saveLevel1a
    disp('Saving Level 1a...')
    calibrationTool = calibrationTool.save_level1a(calibrationTool, logFile,...
        calibratedSpectra, warningLevel0);
end

disp('Warning Level0-1a :')
disp(warningLevel0)

disp('Calibration successful')
calibrationTool.successfulCalibration = true;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%
% Clearing some variables for space
clear rawSpectra; 
clear calibratedSpectra;

% close all force
end
