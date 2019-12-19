function run_retrieval(retrievalTool,dateStr)
% First example for a run function
assert(ischar(dateStr))

% Check here that all required files are filled for retrievalTool !!
% ...

file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', dateStr];
errorLevel0_1a=struct();

try
    [log,rawSpectra]=retrievalTool.read_level0(file);
catch ME
    errorLevel0_1a.readingRawData=ME.identifier;
end

errorLevel0_1a=retrievalTool.check_level0(log,rawSpectra,retrievalTool,errorLevel0_1a);

% 
% switch retrievalTool.instrumentName
%     case 'GROMOS'
%         rawSpectra=retrievalTool.flip_spectra(rawSpectra);
%     case 'SOMORA'
%         % DO NOTHING
% end
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
