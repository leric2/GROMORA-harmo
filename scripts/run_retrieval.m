function run_retrieval(retrievalTool,dateStr)
%%
% First example for a run function
assert(ischar(dateStr),'Please enter the date in the right format')

% Check here that all required fields are filled in retrievalTool !!
retrievalTool_complete(retrievalTool)

file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', dateStr];
assert(exist([file '.txt'],'file') && exist([file '.bin'],'file'),'Files not found')

% Initialize structure containing the error that are non fatal for the
% retrieval
warningLevel0_1a=struct();

% SHOULD WE PUT THE TRY/CATCH INSIDE THE FUNCTIONS ?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading raw data
try
    [log,rawSpectra]=retrievalTool.read_level0(file,retrievalTool);
catch ME
    warningLevel0_1a.readingRawData=ME.identifier;
    error(ME.identifier,'Problem while reading the file')
end

% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible.
if ~strcmp(retrievalTool.instrumentName,'GROMOS') 
    log=retrievalTool.harmonize_log(log);
end

% Quality check of the raw data:
warningLevel0_1a=retrievalTool.check_level0(log,rawSpectra,retrievalTool,warningLevel0_1a);

% Check if the raw data are good
if not(numel(fieldnames(warningLevel0_1a))==0)
    disp('Houston, we have a problem already!')
end

% Reformat the raw spectra from vector to matrix
try
    rawSpectra=retrievalTool.reformat_spectra(rawSpectra,log,retrievalTool);
catch ME
    warningLevel0_1a.reformattingSpectra=ME.identifier;
    error(ME.identifier,'Problem reformatting the spectra')
end

% when needed, flip it !
if retrievalTool.flipped_spectra
    try
        rawSpectra=retrievalTool.flip_spectra(rawSpectra);
    catch ME
        warningLevel0_1a.flippingSpectra=ME.identifier;
        warning(ME.identifier)
    end
end

% Option for plotting spectra (to be improved...)
if retrievalTool.rawSpectraPlot
    try
        retrievalTool.plot_raw_spectra(rawSpectra,0,2e4,20);
    catch ME
        warningLevel0_1a.plottingSpectra=ME.identifier;
        warning(ME.identifier,'s')
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    calibratedSpectra=retrievalTool.calibrate(rawSpectra,log,retrievalTool,80,'time');
catch ME
    warningLevel0_1a.calibrate=ME.identifier;
    warning(ME.identifier,'Problem when doing the calibration')
end

% Option for plotting spectra (to be improved...)
if retrievalTool.calibratedSpectraPlot
    try
        retrievalTool.plot_calibrated_spectra(calibratedSpectra,0,350,25);
    catch ME
        warningLevel0_1a.plottingSpectra=ME.identifier;
        warning(ME.identifier,'Problem Plotting')
    end
end

% Quality check of the calibrated spectra
% Also computing some additional metadata from the log file
try
    calibratedSpectra=retrievalTool.check_calibrated(log,retrievalTool,calibratedSpectra);
catch ME
    warningLevel0_1a.check_calibrated=ME.identifier;
    warning(ME.identifier,'Problem when checking the calibrated spectra')
end

%%
% Saving calibrated spectra (level1a) into NetCDF-4 file
try
    savingLevel0Error=retrievalTool.save_level1a(retrievalTool,calibratedSpectra);
catch ME
    warningLevel0_1a.savingSpectra=ME.identifier;
    warning(ME.identifier,'Problem when saving the calibrated spectra')
end

disp(warningLevel0_1a)
%     
% end
% 
% end
% 
