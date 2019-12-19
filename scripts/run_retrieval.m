function run_retrieval(retrievalTool,dateStr)
% First example for a run function
assert(ischar(dateStr),'Please enter the date in the right format')

% Check here that all required fields are filled in retrievalTool !!
retrievalTool_complete(retrievalTool)

file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', dateStr];
assert(exist([file '.txt'],'file') && exist([file '.bin'],'file'),'Files not found')

% Initialize structure containing the error that are non fatal for the
% retrieval
errorLevel0_1a=struct();

% SHOULD WE PUT THE TRY/CATCH INSIDE THE FUNCTIONS ?
try
    [log,rawSpectra]=retrievalTool.read_level0(file);
catch ME
    errorLevel0_1a.readingRawData=ME.identifier;
end

if ~strcmp(retrievalTool.instrumentName,'GROMOS') 
    log=retrievalTool.harmonize_log(log);
end

errorLevel0_1a=retrievalTool.check_level0(log,rawSpectra,retrievalTool,errorLevel0_1a);

if not(numel(fieldnames(errorLevel0_1a))==0)
    disp('Houston, we have a problem already!')
end

% Reformat the raw spectra from vector to matrix
try
    rawSpectra=retrievalTool.reformat_spectra(rawSpectra,log,retrievalTool);
catch ME
    errorLevel0_1a.reformattingSpectra=ME.identifier;
end

% when needed, flip it !
if retrievalTool.flipped_spectra
    try
        rawSpectra=retrievalTool.flip_spectra(rawSpectra);
    catch ME
        errorLevel0_1a.flippingSpectra=ME.identifier;
    end
end


% 
% 
% try
%     []=retrievalTool.calibrate_data(rawSpectra);
% catch
%     
% end
% 
% try
%     [savingLevel0Error]=retrievalTool.save_calibrated_data(retrievalTool.level1Folder);
% catch
%     
% end
% 
% end
% 
