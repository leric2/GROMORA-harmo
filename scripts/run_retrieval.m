function run_retrieval(retrievalTool)
%==========================================================================
% NAME          | run_retrieval.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      |
%               |
% ABSTRACT      |
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
% Just checking that dateStr is a str...
assert(ischar(retrievalTool.dateStr),'Please enter the date in the right format')

% Check here that all required fields are filled in retrievalTool !!
% Maybe to add later ?
% retrievalTool_complete(retrievalTool)

% Check if both war file exist (does not check their content yet)
assert(exist([retrievalTool.file '.txt'],'file') && exist([retrievalTool.file '.bin'],'file'),'Files not found')
disp(['Starting the calibration process for ' retrievalTool.instrumentName ': ' retrievalTool.dateStr])

if ~retrievalTool.level1aExist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading level0 data...')
% Reading raw data
[logFile,rawSpectra]=retrievalTool.read_level0(retrievalTool);

% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible (different function for each
% instrument and might need date information later ?).
logFile=retrievalTool.harmonize_log(logFile);

% Quality check of the raw data:
if retrievalTool.checkLevel0
    warningLevel0=retrievalTool.check_level0(logFile,rawSpectra,retrievalTool);
else
    warningLevel0='';
end

%% TO CHECK IF RIGHT
% Reformat the raw spectra from vector to matrix
if size(rawSpectra,1)==1
    rawSpectra=retrievalTool.reformat_spectra(rawSpectra,logFile,retrievalTool);
end

% when needed, flip it !
if retrievalTool.flipped_spectra
    rawSpectra=retrievalTool.flip_spectra(rawSpectra);
end

% Option for plotting spectra (to be improved...)
if retrievalTool.rawSpectraPlot
    retrievalTool.plot_raw_spectra(rawSpectra,0,1e9,20);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calibrating...')

% Calibration of the spectra using hot and cold load)
[drift,calibratedSpectra] = retrievalTool.calibrate(rawSpectra,logFile,retrievalTool,retrievalTool.TCold,'standard');

% Quality check of the calibrated spectra
% Also computing some additional metadata from the log file
calibratedSpectra=retrievalTool.check_calibrated(logFile,retrievalTool,calibratedSpectra);

% Option for plotting and saving drift and calibrated spectra
if retrievalTool.calibratedSpectraPlot
    retrievalTool.plot_calibrated_spectra(retrievalTool,drift,calibratedSpectra,50,300,24);
end

%%
% Saving calibrated spectra (level1a) into NetCDF-4 file
disp('Saving Level 1a...')
retrievalTool = retrievalTool.save_level1a(retrievalTool,logFile,calibratedSpectra,warningLevel0);

disp('Warning Level0-1a :')
disp(warningLevel0)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%
% Clearing some variables for space
clear rawSpectra;

clear calibratedSpectra
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1a to level 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We create level1b structure which will contain both the calibrated
% spectra read from the level1a data (level1b.calibratedSpectra and the 
% integrated spectra that will be added later (level1b.integration) 
level1b=struct();

% Defining level1a filename to read (to be adapted for other users)
filename=[retrievalTool.level1Folder retrievalTool.instrumentName '_level1a_' retrievalTool.spectrometer '_' retrievalTool.dateStr '.nc'];
retrievalTool.filenameLevel1a=filename;

if isfield(retrievalTool,'filenameLevel1a') 
    if exist(retrievalTool.filenameLevel1a,'file')
        [level1b.calibratedSpectra,retrievalTool]=retrievalTool.read_level1a(retrievalTool);
    else
        error('The level1a for this file does not exist yet, please calibrate !')
    end
else
    error('No calibration data found for this day')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now because no Payerne dataset
retrievalTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
retrievalTool.get_meteo_data = @(retrievalTool,correctedSpectra) get_meteo_data_unibe(retrievalTool,correctedSpectra);

% Reading meteo data during this day:
level1b.calibratedSpectra=retrievalTool.get_meteo_data(retrievalTool,level1b.calibratedSpectra);

% checking the quality of the channels and flagging the potential bad ones
% (we do not remove any)
level1b.calibratedSpectra=retrievalTool.checking_channel_quality(level1b.calibratedSpectra,retrievalTool,1);

% Compute tropospheric transmittance and correction for every calibrated
% spectra.
level1b.calibratedSpectra=tropospheric_correction_generic(level1b.calibratedSpectra,10.4);

% Integrating the "good spectra" based on tropospheric transmittance and
% calibration flags. --> To improve. Maybe introduce weighted mean of
% spectra based on tropospheric transmittance ?
level1b = retrievalTool.integrate_calibrated_spectra(retrievalTool,level1b);

%% Correction and checks
% Now on the integrated spectra; checking the quality of the channels and 
% flagging the potential bad ones (we do not remove any).
level1b.integration=retrievalTool.checking_channel_quality(level1b.integration,retrievalTool,2);

% Compute tropospheric transmittance and correction for every integrated
% spectra.
level1b.integration=tropospheric_correction_generic(level1b.integration,10.4);

% Performing window correction
level1b=retrievalTool.window_correction(retrievalTool,level1b);

% sideband correction ?
% TODO

% Plotting and saving calibrated and corrected spectra
if retrievalTool.integratedSpectraPlot
    retrievalTool.plot_hourly_spectra(retrievalTool,level1b.integration,50,260)
end

%%
% Saving integrated spectra (level1b) into NetCDF-4 file
disp('Saving Level 1b...')
retrievalTool = retrievalTool.save_level1b(retrievalTool,level1b);

disp('Level 1b saved for : ',retrievalTool.dateStr)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
end
