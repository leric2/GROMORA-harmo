% Main script for launching retrieval
clear; close all; clc;


instrumentName='SOMORA';
dateStart='2019_10_01';
dateEnd='2019_10_02';

dateStr='2019_10_01';

% Import default tools for running a retrieval for a given instrument
% retrievalTool=import_default_retrievalTool(instrumentName);         TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the retrievalTool
retrievalTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};

% Name of the instrument
retrievalTool.instrumentName=instrumentName;

% Valid properties for all instruments
retrievalTool.bytesPerValue=4;

% Path definition (for local computer)
if (instrumentName=='GROMOS')
    retrievalTool.rawFileFolder='/scratch/GROMOS_rawData/2019/';
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/GROMOS/';
    retrievalTool.numberOfChannels=32768;
    retrievalTool.numberOfTippingCurveExpected=48;
    retrievalTool.toleranceTippingCurves=2;
    % Considering the expected number of tipping curve:
    retrievalTool.numberOfCyclesExpected=1500;
    retrievalTool.toleranceNumberCycles=15;
    retrievalTool.tippingSize=27;
    retrievalTool.flipped_spectra=true;
    retrievalTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);
    
elseif (instrumentName=='SOMORA')
    retrievalTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/SOMORA';
    retrievalTool.numberOfChannels=16384;
    retrievalTool.numberOfTippingCurveExpected=48;
    retrievalTool.toleranceTippingCurves=2;
    retrievalTool.tippingSize=5;
    retrievalTool.numberOfCyclesExpected=NaN;
    retrievalTool.toleranceNumberCycles=NaN;
    retrievalTool.flipped_spectra=false;
    % The log needs to be harmonized, for now taking gromos as basis
    retrievalTool.harmonize_log=@(log) harmonize_log_somora(log);
end

% Reading routine to use for the raw data
retrievalTool.read_level0=@(file) read_level0_generic(file);

% Quality check for the raw data
retrievalTool.check_level0= @(log,rawSpectra,retrievalTool,errorLevel0_1a) check_level0_generic(log,rawSpectra,retrievalTool,errorLevel0_1a);

% Reformatting the raw spectra into a matrix (numberOfSpectra x
% numberOfChannels)
retrievalTool.reformat_spectra=@(rawSpectra,log,retrievalTool) reformat_spectra_generic(rawSpectra,log,retrievalTool);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the retrieval with the defined toolchain
run_retrieval(retrievalTool,dateStr)


