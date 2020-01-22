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

file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', retrievalTool.dateStr];
assert(exist([file '.txt'],'file') && exist([file '.bin'],'file'),'Files not found')
disp(['Starting the calibration process for ' retrievalTool.instrumentName ': ' retrievalTool.dateStr])
disp('Reading level0 data...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading raw data
[log,rawSpectra]=retrievalTool.read_level0(file,retrievalTool);

% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible.
log=retrievalTool.harmonize_log(log);

% Quality check of the raw data:
warningLevel0=retrievalTool.check_level0(log,rawSpectra,retrievalTool);

% Reformat the raw spectra from vector to matrix
rawSpectra=retrievalTool.reformat_spectra(rawSpectra,log,retrievalTool);

% when needed, flip it !
if retrievalTool.flipped_spectra
    rawSpectra=retrievalTool.flip_spectra(rawSpectra);
end

% Option for plotting spectra (to be improved...)
if retrievalTool.rawSpectraPlot
    retrievalTool.plot_raw_spectra(rawSpectra,0,2e4,10);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calibrating...')

calibratedSpectra=retrievalTool.calibrate(rawSpectra,log,retrievalTool,80,'standard');

% Quality check of the calibrated spectra
% Also computing some additional metadata from the log file
calibratedSpectra=retrievalTool.check_calibrated(log,retrievalTool,calibratedSpectra);


% Option for plotting and saving spectra (to be improved...)
if retrievalTool.calibratedSpectraPlot
    retrievalTool.plot_calibrated_spectra(retrievalTool,calibratedSpectra,50,350,10);
end

%%
% Clearing some variables for space
clear rawSpectra;

%%
% Saving calibrated spectra (level1a) into NetCDF-4 file
disp('Saving Level 1a...')
retrievalTool.save_level1a(retrievalTool,log,calibratedSpectra,warningLevel0);

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
%
%

% read_level1a()

% plot_hourly_spectra()










