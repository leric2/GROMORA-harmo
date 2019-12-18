% Main script for launching retrieval

instrumentName='GROMOS';
dateStart='2019_10_01';
dateEnd='2019_10_02';

% Import default tools for running a retrieval for a given instrument
% retrievalTool=import_default_retrievalTool(instrumentName);         TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the retrievalTool

% Name of the instrument
retrievalTool.instrumentName=instrumentName;

% Path definition (for local computer)
if (instrumentName=='GROMOS')
    retrievalTool.rawFileFolder='/scratch/GROMOS_RawData/2019/';
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROMOS_SOMORA_HARMO/Level1/GROMOS/';
    
elseif (instrumentName=='SOMORA')
    retrievalTool.rawFileFolder='/scratch/SOMORA_RawData/2019/';
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROMOS_SOMORA_HARMO/Level1/SOMORA';
end

% Reading routine to use for the raw data
retrievalTool.read_level0=@(file) read_level0_generic(file);

% Quality check for the raw data
retrievalTool.check_level0= @(rawSpectra,log) check_level0_generic(rawSpectra,log);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the retrieval with the defined toolchain
run_retrieval(retrievalTool,dateStart)


