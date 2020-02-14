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
% First example for a run function
assert(ischar(retrievalTool.dateStr),'Please enter the date in the right format')

% Check here that all required fields are filled in retrievalTool !!
retrievalTool_complete(retrievalTool)

assert(exist([retrievalTool.file '.txt'],'file') && exist([retrievalTool.file '.bin'],'file'),'Files not found')
disp(['Starting the calibration process for ' retrievalTool.instrumentName ': ' retrievalTool.dateStr])
disp('Reading level0 data...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading raw data
[log,rawSpectra]=retrievalTool.read_level0(retrievalTool);

% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible.
log=retrievalTool.harmonize_log(log);

% Quality check of the raw data:
if retrievalTool.checkLevel0
    warningLevel0=retrievalTool.check_level0(log,rawSpectra,retrievalTool);
else
    warningLevel0='';
end

%% TO CHECK IF RIGHT
% Reformat the raw spectra from vector to matrix
if size(rawSpectra,1)==1
    rawSpectra=retrievalTool.reformat_spectra(rawSpectra,log,retrievalTool);
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

calibratedSpectra = retrievalTool.calibrate(rawSpectra,log,retrievalTool,80,'standard');

% Quality check of the calibrated spectra
% Also computing some additional metadata from the log file
calibratedSpectra=retrievalTool.check_calibrated(log,retrievalTool,calibratedSpectra);

% Option for plotting and saving spectra (to be improved...)
if retrievalTool.calibratedSpectraPlot
    retrievalTool.plot_calibrated_spectra(retrievalTool,calibratedSpectra,50,250,24);
end

%%
% Clearing some variables for space
clear rawSpectra;

%%
% Saving calibrated spectra (level1a) into NetCDF-4 file
disp('Saving Level 1a...')
retrievalTool = retrievalTool.save_level1a(retrievalTool,log,calibratedSpectra,warningLevel0);

disp('Warning Level0-1a :')
disp(warningLevel0)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1a to level 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
%
% 
% if exist(retrievalTool.filenameLevel1a)
%     calibratedSpectra=retrievalTool.read_level1a(retrievalTool);
% end
% 
% % correctedSpectra.date=calibratedSpectra(1).date;
% % 
% 
%calibratedSpectra=retrievalTool.get_meteo_data(calibratedSpectra,retrievalTool);
% 
% 
% % Option for plotting hourly spectra (to be improved...)
% if retrievalTool.hourlyCalibratedSpectraPlot
%     retrievalTool.plot_hourly_spectra(retrievalTool,calibratedSpectra,50,350)
% end
% 
% % checking the quality of the channels and flagging the potential bad ones
% % (we do not remove any)
%calibratedSpectra=retrievalTool.checking_channel_quality(calibratedSpectra,retrievalTool);

% Check calibrated spectra to identify the good ones for integration

% Integrate the "good" spectra --> weighted mean based on tropospheric
% transmittance ?
% retrievalTool.integrate_calibrated_spectra()

% window_correction()

%calibratedSpectra=tropospheric_correction_generic(calibratedSpectra,10.4);

% sideband correction ?

%%
% Saving calibrated spectra (level1a) into NetCDF-4 file
%disp('Saving Level 1b...')
%retrievalTool = retrievalTool.save_level1b(retrievalTool,log,calibratedSpectra,warningLevel0);

%disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')











