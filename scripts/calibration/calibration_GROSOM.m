%==========================================================================
% NAME      | calibration_GROSOM.m
% TYPE      | Script
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 01.2020
%           |
% ABSTRACT  | Main script for launching the calibration of MW
%           | radiometer data. The whole point of this script is to
%           | generate a calibrationTool structure containing every
%           | information needed for the calibration before running it
%           | sequentially. Most of the parameters as imported as
%           | default in the function import_default_calibrationTool()
%           | and the remaining parameters can be defined or modified
%           | here directly. It loops through the dates provided by the
%           | user and launch run_calibration() for each specific
%           | day.
%           |
%           |
% ARGUMENTS | INPUTS: - instrumentName
%           |         - dates: range of date for calibration
%           |
%           | OUTPUTS: Level1a and Level1b, saved from the run_calibration 
%           | function
%           |
% CALLS     | import_default_calibrationTool(instrumentName,dateStr)
%           | run_calibration(calibrationTool)
%           |
%           |
%==========================================================================

clear; close all; clc;

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C'
instrumentName='MIAWARA-C';

% Type of calibration to do: standard of debug
calibrationType='standard';

% Define the dates for the calibration:
dates=datenum('2017_06_08','yyyy_mm_dd'):datenum('2017_06_08','yyyy_mm_dd');
%dates=datenum('2015_09_27','yyyy_mm_dd')

% working directory
root_dir = '/home/franziska/Documents/MW/GROSOM-harmo/';
%cd work_path
addpath(genpath(root_dir))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining all parameters for the calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:numel(dates)
    dateStr=datestr(dates(k),'yyyy_mm_dd');
    
    % Import default tools for running a retrieval for a given instrument
    calibrationTool=import_default_calibrationTool(instrumentName,dateStr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing the calibrationTool for this particular day and instrument:
    % calibrationTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};
    
    % for gaining time.
    calibrationTool.level1aExist=false;
    
    % Time interval for doing the calibration
    calibrationTool.calibrationTime=15;
    
    % Total integration time
    calibrationTool.integrationTime=360; %6h
    
    % Temperature of the cold load
    %calibrationTool.TCold=80;
    
    %
    calibrationTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
    calibrationTool.calendar='standard';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    calibrationTool.hotSpectraNumberOfStdDev=3;
    calibrationTool.coldSpectraNumberOfStdDev=3;

    % Filters for flagging "bad channels"
    % On 10 minutes spectra
    %calibrationTool.filter1.TbMax=300;
    %calibrationTool.filter1.TbMin=20;
    %calibrationTool.filter1.boxCarSize=51;
    %calibrationTool.filter1.boxCarThresh=7;
    
    % On hourly spectra
    %calibrationTool.filter2.TbMax=300;
    %calibrationTool.filter2.TbMin=20;
    %calibrationTool.filter2.boxCarSize=51;
    %calibrationTool.filter2.boxCarThresh=2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug mode and plot options
    
    calibrationTool.calType=calibrationType;
    
    calibrationTool.rawSpectraPlot=true;
    calibrationTool.calibratedSpectraPlot=true;
    calibrationTool.integratedSpectraPlot=true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level0 -> Level1a
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Selecting the functions that will be used for processing this retrieval
    % Reading routine to use for the raw data
%    calibrationTool.read_level0=@(calibrationTool) read_level0_generic(calibrationTool);
    
    % Quality check for the raw data
    calibrationTool.check_level0=@(log,rawSpectra,calibrationTool) check_level0_generic(log,rawSpectra,calibrationTool);
    
    % Reformatting of the raw spectra into a matrix (numberOfSpectra x
    % numberOfChannels)
    calibrationTool.reformat_spectra=@(rawSpectra,log,calibrationTool) reformat_spectra_generic(rawSpectra,log,calibrationTool);
    
    % Plotting some raw spectra:
    calibrationTool.plot_raw_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_raw_spectra_generic(rawSpectra,lowerLim,upperLim,N);
    
    % TODO
    % Find the sky temperature at zenith with a tipping curve
%    calibrationTool.find_T_sky_with_tipping_curve=@(rawSpectra,log,calibrationTool,calType) find_T_sky_with_tipping_curve_generic()
    
    % Function to use for doing the calibration:
%    calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,calType) calibrate_generic(rawSpectra,log,calibrationTool,calType);
    
    % Plot some calibrated spectra:
    calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,rawSpectra,lowerLim,upperLim,N) plot_spectra_generic(calibrationTool,drift,rawSpectra,lowerLim,upperLim,N);
    
    % Function for quality check of the calibrated spectra
%    calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_generic(log,calibrationTool,calibratedSpectra);
    
    % Function saving the calibrated spectra into netCDF file
    calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(calibrationTool,log,calibratedSpectra,warningLevel0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level1a -> Level1b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     % Function reading the daily calibrated spectra from netCDF file
%     calibrationTool.read_level1a = @(calibrationTool) read_level1a_daily(calibrationTool);
%     
%     % Check of the channels quality on the calibrated spectra:
%     calibrationTool.checking_channel_quality= @(calibratedSpectra,calibrationTool,filterN) checking_channel_quality_gromos(calibratedSpectra,calibrationTool,filterN);
%     
%     % Integration of level1a data
%     calibrationTool.integrate_calibrated_spectra= @(calibrationTool,calibratedSpectra) integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra);
%     
%     % Function for plotting the integrated spectra (when hourly)
%     calibrationTool.plot_integrated_spectra = @(calibrationTool,rawSpectra,lowerLim,upperLim) plot_integrated_spectra_generic(calibrationTool,rawSpectra,lowerLim,upperLim);
%     
%     calibrationTool.tropospheric_correction = @(integration,TtropCorr) tropospheric_correction_generic(integration,TtropCorr);
%     
%     % Window correction for the calibrated spectra
%     calibrationTool.window_correction= @(calibrationTool,level1b) window_correction_generic(calibrationTool,level1b);
%     
%     % Function saving the calibrated spectra into netCDF file
%     calibrationTool.save_level1b=@(calibrationTool,level1b) save_level1b_daily(calibrationTool,level1b);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instrument specific parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GROMOS
    if strcmp(instrumentName,'GROMOS')
        % Path definition (for local computer only)
        %calibrationTool.rawFileFolder=['/scratch/GROMOS_rawData/' dateStr(1:4) '/' dateStr(6:7) '/'];
        calibrationTool.rawFileFolder=['/mnt/instrumentdata/gromos/FFTS/' dateStr(1:4) '/'];
        calibrationTool.level1Folder='/home/esauvageat/Documents/GROSOM/Analysis/Level1/GROMOS/';
        calibrationTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
        
        % Some parameters to play with
        calibrationTool.badChannels=[16384 16385];
        
        % Function specific to this instrument
        % meteo Data
        calibrationTool.get_meteo_data = @(calibrationTool,correctedSpectra) get_meteo_data_unibe(calibrationTool,correctedSpectra);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOMORA
    elseif strcmp(instrumentName, 'SOMORA')
        %calibrationTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
        calibrationTool.rawFileFolder=['/home/eric/Documents/PhD/GROSOM/rawData/'];
        %calibrationTool.level1Folder='/home/esauvageat/Documents/GROSOM/Analysis/Level1/SOMORA/';
        calibrationTool.level1Folder='/home/eric/Documents/PhD/GROSOM/Level1/';
        calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
        % TOCHANGE
        %calibrationTool.meteoFolder='/home/esauvageat/Documents/GROSOM/Analysis/MeteoFile/METEO_DATA/';
        calibrationTool.meteoFolder='/home/eric/Documents/PhD/METEO_DATA/';
        
        
        calibrationTool.badChannels=1:104;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function specific to this instrument
        % meteo Data
        calibrationTool.get_meteo_data = @(calibrationTool,correctedSpectra) get_meteo_data_payerne(calibrationTool,correctedSpectra);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MOPI
    elseif strcmp(instrumentName,'mopi5')
        % FOR MOPI:
        % Everything stored into "import_default_calibrationTool"
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MIAWARA-C
    elseif strcmp(instrumentName,'MIAWARA-C')
        % FOR MIAWARA-C:
        % Everything stored into "import_default_calibrationTool"
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launching the calibration process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now, we keep it dirty for separating between GROSOM and MOPI
    if calibrationTool.numberOfSpectrometer==1
        try
            % if commented, nothing happens --> developping purposes
            % run_calibration(calibrationTool)
        end
    else
        for m=1:3
            model=[1 3 4];
            fftsMopi=model(m);
            calibrationTool.ffts_model=fftsMopi;
            S  = {'USRP-A', 'USRP-B','U5303', 'AC240'};
            calibrationTool.spectrometer=S{calibrationTool.ffts_model};
            try  
                % Running the retrieval with the defined toolchain
                %run_calibration(calibrationTool)
            end
        end       
    end
end

