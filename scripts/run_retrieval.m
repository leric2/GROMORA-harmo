function run_retrieval(retrievalTool,dateStr)
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
try
    [log,rawSpectra]=retrievalTool.read_level0(file,retrievalTool);
catch ME
    warningLevel0_1a.readingRawData=ME.identifier;
    error(ME.identifier,'Problem while reading the file')
end

if ~strcmp(retrievalTool.instrumentName,'GROMOS') 
    log=retrievalTool.harmonize_log(log);
end

warningLevel0_1a=retrievalTool.check_level0(log,rawSpectra,retrievalTool,warningLevel0_1a);

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
        retrievalTool.plot_spectra(rawSpectra,0,2e4,20);
    catch ME
        warningLevel0_1a.plottingSpectra=ME.identifier;
        warning(ME.identifier,'s')
    end
end

disp(warningLevel0_1a)

try
    calibratedSpectra=retrievalTool.calibrate(rawSpectra,log,retrievalTool,80,'time');
catch
    warningLevel0_1a.calibrate=ME.identifier;
    warning(ME.identifier)
end

% Option for plotting spectra (to be improved...)
if retrievalTool.calibratedSpectraPlot
    try
        retrievalTool.plot_spectra(calibratedSpectra.Tb,100,350,25);
    catch ME
        warningLevel0_1a.plottingSpectra=ME.identifier;
        warning(ME.identifier,'Problem Plotting')
    end
end

% 
% try
%     [savingLevel0Error]=retrievalTool.save_calibrated_data(retrievalTool.level1Folder);
% catch
%     
% end
% 
% end
% 
