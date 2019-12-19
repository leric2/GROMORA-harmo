% Main script for launching retrieval
clear; close all; clc;


instrumentName='GROMOS';
dateStart='2019_10_01';
dateEnd='2019_10_02';

dateStr='2019_10_01';

% Import default tools for running a retrieval for a given instrument
% retrievalTool=import_default_retrievalTool(instrumentName);         TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the retrievalTool

% Name of the instrument
retrievalTool.instrumentName=instrumentName;

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
    
elseif (instrumentName=='SOMORA')
    retrievalTool.rawFileFolder='/scratch/SOMORA_rawData/2019/';
    retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/SOMORA';
    retrievalTool.numberOfChannels=16384;
    retrievalTool.numberOfTippingCurveExpected=48;
    retrievalTool.tippingSize=5;
end

% Reading routine to use for the raw data
retrievalTool.read_level0=@(file) read_level0_generic(file);

% Quality check for the raw data
retrievalTool.check_level0= @(log,rawSpectra,retrievalTool,errorLevel0_1a) check_level0_generic(log,rawSpectra,retrievalTool,errorLevel0_1a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the retrieval with the defined toolchain
run_retrieval(retrievalTool,dateStart)


