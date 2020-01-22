%==========================================================================
% NAME          | retrieval.m
% TYPE          | Main script for launching the retrieval
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 12.2019
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
clear; close all; clc;

instrumentName='SOMORA';

% Define the dates where we want to launch a retrieval:
dates=datenum('2019_08_01','yyyy_mm_dd'):datenum('2019_08_05','yyyy_mm_dd');

for k = 1:numel(dates)
    dateStr=datestr(dates(k),'yyyy_mm_dd');
    
    % Import default tools for running a retrieval for a given instrument
    retrievalTool=import_default_retrievalTool(instrumentName);
    
    retrievalTool.dateStr=dateStr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing the retrievalTool for this particular retrieval:
    retrievalTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};
    
    retrievalTool.calibrationTime=10;
    
    % Path definition (for local computer)
    if (instrumentName=='GROMOS')
        retrievalTool.rawFileFolder=['/scratch/GROMOS_rawData/' dateStr(1:4) '/' dateStr(6:7) '/'];
        %retrievalTool.rawFileFolder=['/mnt/instrumentdata/gromos/FFTS/' dateStr(1:4) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/GROMOS/';
        
    elseif (instrumentName=='SOMORA')
        retrievalTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/SOMORA/';
    end
    
    retrievalTool.saveAllCycles=0;
    
    retrievalTool.hotSpectraNumberOfStdDev=3;
    retrievalTool.coldSpectraNumberOfStdDev=3;
    
    retrievalTool.rawSpectraPlot=false;
    retrievalTool.calibratedSpectraPlot=true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Selecting the functions that will be used for processing this retrieval
    %
    % Reading routine to use for the raw data
    retrievalTool.read_level0=@(file,retrievalTool) read_level0_generic(file,retrievalTool);
    
    % Quality check for the raw data
    retrievalTool.check_level0=@(log,rawSpectra,retrievalTool,errorLevel0_1a) check_level0_generic(log,rawSpectra,retrievalTool,errorLevel0_1a);
    
    % Reformatting of the raw spectra into a matrix (numberOfSpectra x
    % numberOfChannels)
    retrievalTool.reformat_spectra=@(rawSpectra,log,retrievalTool) reformat_spectra_generic(rawSpectra,log,retrievalTool);
    
    % Plotting some raw spectra:
    retrievalTool.plot_raw_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_raw_spectra_generic(rawSpectra,lowerLim,upperLim,N);
    
    retrievalTool.plot_calibrated_spectra=@(retrievalTool,rawSpectra,lowerLim,upperLim,N) plot_spectra_generic(retrievalTool,rawSpectra,lowerLim,upperLim,N);
    
    retrievalTool.calibrate=@(rawSpectra,log,retrievalTool,TCold,calType) calibrate_generic(rawSpectra,log,retrievalTool,TCold,calType);
    
    retrievalTool.check_calibrated=@(log,retrievalTool,calibratedSpectra) check_calibrated_generic(log,retrievalTool,calibratedSpectra);
    
    retrievalTool.save_level1a=@(retrievalTool,log,calibratedSpectra) save_level1a_daily(retrievalTool,log,calibratedSpectra);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Running the retrieval with the defined toolchain
    run_retrieval(retrievalTool)
end

