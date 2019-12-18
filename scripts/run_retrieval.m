function run_retrieval(retrievalTool,date)
% First example for a run function
assert(ischar(date))

file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09', date];

try
    [log,rawSpectra,readingLevel0Error]=retrievalTool.read_level0(file);
catch
    
end
% 
% try
%     [qualityLevel0]=retrievalTool.check_level0(log,rawSpectra);
% catch
%     
% end
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
