%==========================================================================
% NAME          | retrieval.m
% TYPE          | Script
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Main script for launching the calibration of MW
%               | radiometer data. The whole point of this script is to
%               | generate a retrievalTool structure containing every
%               | information needed for the calibration before running it
%               | sequentially. 
%               | 
%               |
%               |
% ARGUMENTS     | 
%               |
%               | 
%              11
% CALLS         | run_retrieval.m
%               |
%               |
%               |

%==========================================================================
clear; close all; clc;

% 'GROMOS' // 'SOMORA' // 'mopi5'
instrumentName='mopi5';

% Define the dates where we want to launch a retrieval:
dates=datenum('2019_06_16','yyyy_mm_dd'):datenum('2019_06_16','yyyy_mm_dd');

for k = 1:numel(dates)
    dateStr=datestr(dates(k),'yyyy_mm_dd');
    
    % Import default tools for running a retrieval for a given instrument
    retrievalTool=import_default_retrievalTool(instrumentName);
    
    retrievalTool.dateStr=dateStr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing the retrievalTool for this particular retrieval:
    retrievalTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};
    
    retrievalTool.calibrationTime=10;
    
    retrievalTool.observationFreq=1.4217504e11;
    
    retrievalTool.saveAllCycles=1;
    
    retrievalTool.hotSpectraNumberOfStdDev=3;
    retrievalTool.coldSpectraNumberOfStdDev=3;
    
    retrievalTool.rawSpectraPlot=false;
    retrievalTool.calibratedSpectraPlot=true;
    retrievalTool.hourlyCalibratedSpectraPlot=true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Selecting the functions that will be used for processing this retrieval
    %
    % Reading routine to use for the raw data
    retrievalTool.read_level0=@(retrievalTool) read_level0_generic(retrievalTool);
    
    % Quality check for the raw data
    retrievalTool.check_level0=@(log,rawSpectra,retrievalTool) check_level0_generic(log,rawSpectra,retrievalTool);
    
    % Reformatting of the raw spectra into a matrix (numberOfSpectra x
    % numberOfChannels)
    retrievalTool.reformat_spectra=@(rawSpectra,log,retrievalTool) reformat_spectra_generic(rawSpectra,log,retrievalTool);
    
    % Plotting some raw spectra:
    retrievalTool.plot_raw_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_raw_spectra_generic(rawSpectra,lowerLim,upperLim,N);
    
    % Function to use for doing the calibration:
    retrievalTool.calibrate=@(rawSpectra,log,retrievalTool,TCold,calType) calibrate_generic(rawSpectra,log,retrievalTool,TCold,calType);
    
    % Plot some calibrated spectra:
    retrievalTool.plot_calibrated_spectra=@(retrievalTool,rawSpectra,lowerLim,upperLim,N) plot_spectra_generic(retrievalTool,rawSpectra,lowerLim,upperLim,N);
    
    % Function for quality check of the calibrated spectra
    retrievalTool.check_calibrated=@(log,retrievalTool,calibratedSpectra) check_calibrated_generic(log,retrievalTool,calibratedSpectra);
    
    % Function saving the calibrated spectra into netCDF file
    retrievalTool.save_level1a=@(retrievalTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(retrievalTool,log,calibratedSpectra,warningLevel0);
    
    % Function reading the daily calibrated spectra from netCDF file
    retrievalTool.read_level1a = @(retrievalTool) read_level1a_daily(retrievalTool);
    
    %
    retrievalTool.plot_hourly_spectra = @(retrievalTool,rawSpectra,lowerLim,upperLim) plot_hourly_spectra_generic(retrievalTool,rawSpectra,lowerLim,upperLim);
    
    % TO CHANGE FOR SOMORA
    retrievalTool.get_meteo_data = @(retrievalTool,correctedSpectra) get_meteo_data_unibe(retrievalTool,correctedSpectra);
    
    retrievalTool.checking_channel_quality= @(calibratedSpectra,retrievalTool) checking_channel_quality_gromos(calibratedSpectra,retrievalTool);
    
    % Path definition (for local computer)
    if (string(instrumentName)=='GROMOS')
        %retrievalTool.rawFileFolder=['/scratch/GROMOS_rawData/' dateStr(1:4) '/' dateStr(6:7) '/'];
        retrievalTool.rawFileFolder=['/mtn/datalakeMW/instrumentdata/gromos/FFTS/' dateStr(1:4) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/GROMOS/';
        retrievalTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        retrievalTool.file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', retrievalTool.dateStr];
        
        retrievalTool.fLO1=1.45875e11;
        retrievalTool.fLO2=3.6e9;
        
        % This one should correspond to the DC channel
        retrievalTool.LOFreqTot=retrievalTool.fLO1-retrievalTool.fLO2;
        retrievalTool.DCChannel=16384; %=Nchannel/2 ??
        
    elseif (string(instrumentName)=='SOMORA')
        retrievalTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/SOMORA/';
        retrievalTool.file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', retrievalTool.dateStr];
        % TOCHANGE
        retrievalTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        
        retrievalTool.fLO1=1.49275e11;
        retrievalTool.fLO2=5.6e9;
        retrievalTool.fLO3=2e9;
        
        % This one should correspond to the DC channel
        retrievalTool.LOFreqTot=retrievalTool.fLO1-retrievalTool.fLO2-retrievalTool.fLO3;
        retrievalTool.DCChannel=1; %=Nchannel/2 ??
        
    elseif (string(instrumentName)=='mopi5')
        retrievalTool.rawFileFolder=['/mtn/datalakeMW/instrumentdata/mopi5/' dateStr(1:4) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/MOPI5/Level1/';
        
        retrievalTool.file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'_', retrievalTool.dateStr(1:4) retrievalTool.dateStr(6:7) retrievalTool.dateStr(9:10)];
        
        retrievalTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        retrievalTool.observationFreq=110;
        
        %retrievalTool.fLO1=1.49275e11;
        %retrievalTool.fLO2=5.6e9;
        %retrievalTool.fLO3=2e9;
        
        % This one should correspond to the DC channel
        %retrievalTool.LOFreqTot=retrievalTool.fLO1-retrievalTool.fLO2-retrievalTool.fLO3;
        %retrievalTool.DCChannel=1; %=Nchannel/2 ??
        retrievalTool.ffts_model=3;
        retrievalTool.read_level0=@(retrievalTool) mopi5_read(retrievalTool); 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Running the retrieval with the defined toolchain
    %run_retrieval(retrievalTool)
end

