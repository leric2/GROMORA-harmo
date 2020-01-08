% Main script for launching retrieval
clear; close all; clc;


instrumentName='SOMORA';
dateStart='2019_10_01';
dateEnd='2019_10_02';

dateStr='2019_10_03';

% Import default tools for running a retrieval for a given instrument
retrievalTool=import_default_retrievalTool(instrumentName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the retrievalTool for this particular retrieval:
retrievalTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};

% Valid properties for all instruments
retrievalTool.rawSpectraPlot=false;
retrievalTool.calibrationTime=30;
retrievalTool.calibratedSpectraPlot=true;

% Path definition (for local computer)
if (instrumentName=='GROMOS')
    retrievalTool.rawFileFolder='/scratch/GROMOS_rawData/2019/';
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/GROMOS/';
    
elseif (instrumentName=='SOMORA')
    retrievalTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/SOMORA/';
end

% Reading routine to use for the raw data
retrievalTool.read_level0=@(file,retrievalTool) read_level0_generic(file,retrievalTool);

% Quality check for the raw data
retrievalTool.check_level0=@(log,rawSpectra,retrievalTool,errorLevel0_1a) check_level0_generic(log,rawSpectra,retrievalTool,errorLevel0_1a);

% Reformatting the raw spectra into a matrix (numberOfSpectra x
% numberOfChannels)
retrievalTool.reformat_spectra=@(rawSpectra,log,retrievalTool) reformat_spectra_generic(rawSpectra,log,retrievalTool);

retrievalTool.plot_raw_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_raw_spectra_generic(rawSpectra,lowerLim,upperLim,N);

retrievalTool.plot_calibrated_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_spectra_generic(rawSpectra,lowerLim,upperLim,N);

retrievalTool.calibrate=@(rawSpectra,log,retrievalTool,TCold,calType) calibrate_generic(rawSpectra,log,retrievalTool,TCold,calType);

retrievalTool.check_calibrated=@(log,retrievalTool,calibratedSpectra) check_calibrated_generic(log,retrievalTool,calibratedSpectra);

retrievalTool.save_level1a=@(retrievalTool,calibratedSpectra,dateStr) save_level1a_generic(retrievalTool,calibratedSpectra,dateStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the retrieval with the defined toolchain
run_retrieval(retrievalTool,dateStr)


