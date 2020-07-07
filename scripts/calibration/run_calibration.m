function run_calibration(calibrationTool)
%==========================================================================
% NAME      | run_calibration.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 01.2020
%           |
% ABSTRACT  | The main function executing the calibration for the
%           | instrument defined in calibrationTool. Some parts and
%           | functions are dependent on the instrument that we want to
%           | calibrate. 
%           | 
%           |
% ARGUMENTS | INPUTS: - calibrationTool: structure containing all
%           | information about the calibration we want to perform.
%           | Documentation about this structure can be found in external
%           | document.
%           |
%           |
%           | OUTPUTS: - Level1a
%           |          - Level1b 
%           |
%           |
% CALLS     | Some depends on instruments, all are stored in calibrationTool:
%           | %%%%%%%%%%%%%%%%%%%%% Level0 -> Level1a
%           | read_level0(calibrationTool)
%           | harmonize_log(logFile)
%           | check_level0
%           | reformat_spectra
%           | flip_spectra
%           | plot_raw_spectra
%           | calibrate
%           | check_calibrated
%           | plot_calibrated_spectra
%           | save_level1a
%           | %%%%%%%%%%%%%%%%%%%%% Level1a -> Level1b
%           | read_level1a
%           | get_meteo_data
%           | checking_channel_quality
%           | tropospheric_correction_generic
%           | window_correction
%           | plot_integrated_spectra
%           | save_level1b
%           |
%==========================================================================
% Just checking that dateStr is a str...
assert(ischar(calibrationTool.dateStr),'Please enter the date in the right format')

% Check here that all required fields are filled in calibrationTool !!
% Check if needed ?
% calibrationTool_complete(calibrationTool)

% Check if both raw file exist (does not check their content yet)
assert(exist([calibrationTool.file '.txt'],'file') && exist([calibrationTool.file '.bin'],'file'),'Files not found')

% Start calibration
disp(['Starting the calibration process for ' calibrationTool.instrumentName ': ' calibrationTool.dateStr])
%%
if ~calibrationTool.level1aExist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading level0 data...')

% Reading raw data
[logFile,rawSpectra] = calibrationTool.read_level0(calibrationTool);

% Check the temporal content of the file
% TODO 

% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible (different function for each
% instrument and might need date information later ?).
logFile = calibrationTool.harmonize_log(logFile);

% Quality check of the raw data:
if calibrationTool.checkLevel0
    warningLevel0 = calibrationTool.check_level0(logFile,rawSpectra,calibrationTool);
    
    %==============> write overflow spectra
else
    warningLevel0 = '';
end

%% TO CHECK IF RIGHT
% Reformat the raw spectra from vector to matrix
if size(rawSpectra,1) == 1
    rawSpectra = calibrationTool.reformat_spectra(rawSpectra,logFile,calibrationTool);
end

% when needed, flip it !
if calibrationTool.flipped_spectra
    rawSpectra = calibrationTool.flip_spectra(rawSpectra);
end

% TODO
% remove channels which are known to be bad (2pol)
% do not delete channels yet, -> set high error during retrieval

% TODO
% check_level0


% Option for plotting spectra (to be improved...)
if calibrationTool.rawSpectraPlot
    calibrationTool.plot_raw_spectra(rawSpectra,0,1e4,20);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_meteo_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%logFile.meteo = calibrationTool.get_meteo_data(calibrationTool);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tipping Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%logFile.TC = calibrationTool.run_tipping_curve(rawSpectra, logFile, calibrationTool);


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
[drift,calibratedSpectra] = calibrationTool.calibrate(rawSpectra,logFile,calibrationTool,calibrationTool.calType);

% Quality check of the calibrated spectra
% Also computing some additional metadata from the log file and storing
% everything in calibrated spectra
calibratedSpectra = calibrationTool.check_calibrated(logFile,calibrationTool,calibratedSpectra);

% Option for plotting and saving drift and calibrated spectra
if calibrationTool.calibratedSpectraPlot
    try
        calibrationTool.plot_calibrated_spectra(calibrationTool,drift,calibratedSpectra,50,300,24);
    catch ME
        warning(ME.identifier,'problem with the plotting:');
        disp(ME.message)
    end
end

% Saving calibrated spectra (level1a) into NetCDF-4 file
disp('Saving Level 1a...')
calibrationTool = calibrationTool.save_level1a(calibrationTool,logFile,calibratedSpectra,warningLevel0);

disp('Warning Level0-1a :')
disp(warningLevel0)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%
% Clearing some variables for space
clear rawSpectra; 
clear calibratedSpectra;

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1a to level 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We create level1b structure which will contain both the calibrated
% spectra read from the level1a data (level1b.calibratedSpectra and the 
% integrated spectra that will be added later (level1b.integration) 
level1b = struct();

% Defining level1a filename to read (to be adapted for other users)
filename = [calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];
calibrationTool.filenameLevel1a = filename;

if isfield(calibrationTool,'filenameLevel1a') 
    if exist(calibrationTool.filenameLevel1a,'file')
        [level1b.calibratedSpectra,calibrationTool] = calibrationTool.read_level1a(calibrationTool);
    else
        error('The level1a for this file does not exist yet, please calibrate !')
    end
else
    error('No calibration data found for this day')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now because no Payerne dataset
% calibrationTool.meteoFolder = '/mnt/instrumentdata/meteo/exwi/meteo/';
% calibrationTool.get_meteo_data  =  @(calibrationTool,correctedSpectra) get_meteo_data_unibe(calibrationTool,correctedSpectra);

% Reading meteo data during this day:
level1b.calibratedSpectra = calibrationTool.get_meteo_data(calibrationTool,level1b.calibratedSpectra);

% checking the quality of the channels and flagging the potential bad ones
% (we do not remove any)
level1b.calibratedSpectra = calibrationTool.checking_channel_quality(level1b.calibratedSpectra,calibrationTool,1);

% Compute tropospheric transmittance and correction for every calibrated
% spectra.
level1b.calibratedSpectra = calibrationTool.tropospheric_correction(level1b.calibratedSpectra,10.4);

% Integrating the "good spectra" based on tropospheric transmittance and
% calibration flags. --> To improve. Maybe introduce weighted mean of
% spectra based on tropospheric transmittance ?
level1b = calibrationTool.integrate_calibrated_spectra(calibrationTool,level1b);

%% Correction and checks
% Now on the integrated spectra; checking the quality of the channels and 
% flagging the potential bad ones (we do not remove any).
level1b.integration = calibrationTool.checking_channel_quality(level1b.integration,calibrationTool,2);

% Performing window correction
level1b = calibrationTool.window_correction(calibrationTool,level1b);

% Compute tropospheric transmittance and correction for every integrated
% spectra.
level1b.integration = calibrationTool.tropospheric_correction(level1b.integration,10.4);

% sideband correction ?
% TODO

% Plotting and saving calibrated and corrected spectra
if calibrationTool.integratedSpectraPlot
    calibrationTool.plot_integrated_spectra(calibrationTool,level1b.integration,50,260)
end

%%
% Saving integrated spectra (level1b) into NetCDF-4 file
disp('Saving Level 1b...')
calibrationTool  =  calibrationTool.save_level1b(calibrationTool,level1b);

end
