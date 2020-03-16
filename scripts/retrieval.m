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
% ARGUMENTS     | NONE
%               |
%               | 
%               |
% CALLS         | run_retrieval(retrievalTool);
%               |
%               |
%               |

%==========================================================================
clear; close all; clc;

% 'GROMOS' // 'SOMORA' // 'mopi5'
instrumentName='SOMORA';

% Define the dates for the calibration:
dates=datenum('2019_04_16','yyyy_mm_dd'):datenum('2019_04_16','yyyy_mm_dd');
k=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:numel(dates)
    dateStr=datestr(dates(k),'yyyy_mm_dd');
    
    % Import default tools for running a retrieval for a given instrument
    retrievalTool=import_default_retrievalTool(instrumentName,dateStr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing the retrievalTool for this particular day and instrument:
    % retrievalTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};
    
    % for gaining time..
    retrievalTool.level1aExist=true;
    
    % Time interval for doing the calibration
    retrievalTool.calibrationTime=10;
    
    % Total integration time
    retrievalTool.integrationTime=60;
    
    % Temperature of the cold load
    retrievalTool.TCold=80;
    
    %
    retrievalTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
    retrievalTool.calendar='standard';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    retrievalTool.hotSpectraNumberOfStdDev=3;
    retrievalTool.coldSpectraNumberOfStdDev=3;
    
    
    % Filters for flagging "bad channels"
    % On 10 minutes spectra
    retrievalTool.filter1.TbMax=300;
    retrievalTool.filter1.TbMin=20;
    retrievalTool.filter1.boxCarSize=51;
    retrievalTool.filter1.boxCarThresh=7;
    
    % On hourly spectra
    retrievalTool.filter2.TbMax=300;
    retrievalTool.filter2.TbMin=20;
    retrievalTool.filter2.boxCarSize=51;
    retrievalTool.filter2.boxCarThresh=2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug mode and plot options
    retrievalTool.saveAllCycles=0;
    
    retrievalTool.rawSpectraPlot=false;
    retrievalTool.calibratedSpectraPlot=true;
    retrievalTool.integratedSpectraPlot=true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level0 -> Level1a
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Selecting the functions that will be used for processing this retrieval
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
    retrievalTool.plot_calibrated_spectra=@(retrievalTool,drift,rawSpectra,lowerLim,upperLim,N) plot_spectra_generic(retrievalTool,drift,rawSpectra,lowerLim,upperLim,N);
    
    % Function for quality check of the calibrated spectra
    retrievalTool.check_calibrated=@(log,retrievalTool,calibratedSpectra) check_calibrated_generic(log,retrievalTool,calibratedSpectra);
    
    % Function saving the calibrated spectra into netCDF file
    retrievalTool.save_level1a=@(retrievalTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(retrievalTool,log,calibratedSpectra,warningLevel0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level1a -> Level1b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Function reading the daily calibrated spectra from netCDF file
    retrievalTool.read_level1a = @(retrievalTool) read_level1a_daily(retrievalTool);
    
    %
    retrievalTool.plot_hourly_spectra = @(retrievalTool,rawSpectra,lowerLim,upperLim) plot_hourly_spectra_generic(retrievalTool,rawSpectra,lowerLim,upperLim);
    
    % Check of the channels quality on the calibrated spectra:
    retrievalTool.checking_channel_quality= @(calibratedSpectra,retrievalTool,filterN) checking_channel_quality_gromos(calibratedSpectra,retrievalTool,filterN);
    
    % Integration of level1a data
    retrievalTool.integrate_calibrated_spectra= @(retrievalTool,calibratedSpectra) integrate_calibrated_spectra_generic(retrievalTool,calibratedSpectra);
    
    % Window correction for the calibrated spectra
    retrievalTool.window_correction= @(retrievalTool,level1b) window_correction_generic(retrievalTool,level1b);
    
    % Function saving the calibrated spectra into netCDF file
    retrievalTool.save_level1b=@(retrievalTool,level1b) save_level1b_daily(retrievalTool,level1b);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GROMOS
    if (string(instrumentName)=='GROMOS')
        % Path definition (for local computer only)
        %retrievalTool.rawFileFolder=['/scratch/GROMOS_rawData/' dateStr(1:4) '/' dateStr(6:7) '/'];
        retrievalTool.rawFileFolder=['/mnt/instrumentdata/gromos/FFTS/' dateStr(1:4) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/GROMOS/';
        retrievalTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        retrievalTool.file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', retrievalTool.dateStr];
        
        
        % Some parameters to play with
        retrievalTool.badChannels=[16384 16385];
        
        % Function specific to this instrument
        % meteo Data
        retrievalTool.get_meteo_data = @(retrievalTool,correctedSpectra) get_meteo_data_unibe(retrievalTool,correctedSpectra);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOMORA
    elseif (string(instrumentName)=='SOMORA')
        retrievalTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/GROSOM/Level1/SOMORA/';
        retrievalTool.file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'09_', retrievalTool.dateStr];
        % TOCHANGE
        retrievalTool.meteoFolder='/home/esauvageat/Documents/GROSOM/MeteoFile/';
        
        
        retrievalTool.badChannels=1:104;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function specific to this instrument
        % meteo Data
        retrievalTool.get_meteo_data = @(retrievalTool,correctedSpectra) get_meteo_data_payerne(retrievalTool,correctedSpectra);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MOPI
    elseif (string(instrumentName)=='mopi5')
        % FOR MOPI:
        % Everything stored into "import_default_retrievalTool"
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launching the processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if retrievalTool.numberOfSpectrometer==1
        try
            %run_retrieval(retrievalTool)
        end
    else
        for m=1:3
            model=[1 3 4];
            fftsMopi=model(m);
            retrievalTool.ffts_model=fftsMopi;
            S  = {'USRP-A', 'USRP-B','U5303', 'AC240'};
            retrievalTool.spectrometer=S{retrievalTool.ffts_model};
            try  
                % Running the retrieval with the defined toolchain
                %run_retrieval(retrievalTool)
            end
        end       
    end
end

