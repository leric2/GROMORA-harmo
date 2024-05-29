function calibration_GROMORA_reprocessing(dateStr)
% NAME      | calibration_GROMORA_reprocessing.m
% TYPE      | Script
% AUTHOR(S) | Alistair Bell
% CREATION  | 05.2024
%           |
% ABSTRACT  | This script, based on Eric Sauvageat's 'calibration_GROMORA',
%           | extends its functionality to accept a string date argument
%           | ('yyyy_mm_dd'). It performs calibration of gromos
%           | data, giving outputs in Level1a and Level1b formats.
%           | This adaptation allows for targeted reprocessing of data for 
%           | specific dates.
%
% USAGE     | Call with a date stringe, e.g.: 'calibration_GROMORA_reprocessing('2024_05_29')'.


clear all; close all; clc; clear functions; %clear mex;

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C' // 'WIRAC' //
instrumentName='GROMOS';

% Type of calibration to do: standard or debug
calibrationType='standard';

calibrate = true;
integrate = true;
readLabviewLog = true;

% GROMOS from 10.03.2010 only (after change in SW, see logfile), meteo from
customDate = datenum(dateStr, 'yyyy_mm_dd');
dates = customDate : customDate;

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C' // 'WIRAC' //
instrumentName='GROMOS';

% Type of calibration to do: standard or debug
calibrationType='standard';

calibrate = true;
integrate = true;
readLabviewLog = true;


if (strcmp(instrumentName,'GROMOS') | strcmp(instrumentName,'SOMORA')) & readLabviewLog
    labviewLogFolder = '/storage/atmosphere/instruments/gromos/level1/GROMORA/InputsCalibration/';
    labviewLog = read_labview_log_generic(instrumentName, labviewLogFolder);
else
    labviewLog = struct();
end

% % working directory
% %root_dir = '/home/franziska/Documents/MW/GROSOM-harmo/';
% root_dir = pwd;
% %cd work_path
% addpath(genpath(root_dir))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining all parameters for the calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for d = 1:numel(dates)
    dateStr=datestr(dates(d),'yyyy_mm_dd');
    
    % Import default tools for running a calibration or integration for a 
    % given instrument
    calibrationTool=import_default_calibrationTool(dateStr);
    
    calibrationTool.instrumentName=instrumentName;
    calibrationTool.calibrationVersion = '2.0';
    
    calibrationTool.extraName = '';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing the calibrationTool for this particular day and instrument:
    % calibrationTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};
    
    calibrationTool.meanDatetimeUnit='days since 2000-01-01 00:00:00';
    calibrationTool.referenceTime = datenum(2000,1,1,0,0,0);
    calibrationTool.calendar='standard';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    calibrationTool.hotSpectraNumberOfStdDev=3;
    calibrationTool.coldSpectraNumberOfStdDev=3;
    calibrationTool.skySpectraNumberOfStdDev=6;
    
    calibrationTool.labviewLog = labviewLog;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug mode and plot options
    calibrationTool.calType=calibrationType;
    
    % Option for outlier detection during calibration (only with standard
    % calibration type):
    % 'standard', 'noFFT', 'none'
    calibrationTool.outlierDectectionType = 'standard';
    
    calibrationTool.rawSpectraPlot = false;
    calibrationTool.calibratedSpectraPlot = true;
    calibrationTool.calibratedSpectraSpectralPlot = true;
    calibrationTool.saveLevel1a = true;
    calibrationTool.integratedSpectraPlot = true;
    calibrationTool.saveLevel1b = true;
    
    calibrationTool.savePlanckIntensity = true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instrument specific parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each instruments, we define some key parameters here and we
    % import the rest of the parameters and function from their specific
    % import calibrationTool routines. 
    
    if strcmp(instrumentName,'GROMOS')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        
        % Filtering options of level 1a
        calibrationTool.filterByTransmittance = true;
        calibrationTool.transmittanceThreshold = 0.05;
        calibrationTool.filterByFlags = true;
        
        %%% Flags level 1b:
        % Minimum number of averaged spectra needed for the level 1b
        calibrationTool.minNumberOfAvgSpectra = 3;
        % transmittance threshold for flagging level 1b:
        calibrationTool.troposphericTransmittanceFlag = 0.15;

        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % Import specific parameter and functions for this instrument
        calibrationTool = import_GROMOS_calibrationTool(calibrationTool);
        
        % an extra folder where we copy missing anetz data (from STARTWAVE)
        calibrationTool.meteoAnetzExtraFolder = '/storage/atmosphere/MeteoSchweiz/extra/';

    elseif strcmp(instrumentName, 'SOMORA')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        
        % Filtering options of level 1a
        calibrationTool.filterByTransmittance = true;
        calibrationTool.transmittanceThreshold = 0.05;
        calibrationTool.filterByFlags = true;
        
        %%% Flags level 1b:
        % Minimum number of averaged spectra needed for the level 1b
        calibrationTool.minNumberOfAvgSpectra = 3;
        % transmittance threshold for flagging level 1b:
        calibrationTool.troposphericTransmittanceFlag = 0.15;

        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % Import specific parameter and functions for this instrument
        calibrationTool = import_SOMORA_calibrationTool(calibrationTool);

        % an extra folder where we copy missing anetz data (from STARTWAVE)
        calibrationTool.meteoAnetzExtraFolder = '/storage/atmosphere/MeteoSchweiz/extra/';

    elseif strcmp(instrumentName,'mopi5')
        % FOR MOPI:
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=6*60;
        
        calibrationTool.filterByTransmittance = false; % Best to keep false for MOPI5 studies
        calibrationTool.transmittanceThreshold = 0.2;
        calibrationTool.filterByFlags = false;
        
        calibrationTool.minNumberOfAvgSpectra = 6;
        calibrationTool.troposphericTransmittanceFlag = 0.2;
        
        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % the number of the spectrometer models we are interested in
        % see order in calibrationTool.spectrometerTypes
        modelFFTS=[1 3 4];
        
        calibrationTool = import_MOPI5_calibrationTool(calibrationTool);
    elseif strcmp(instrumentName,'MIAWARA-C')
        % FOR MIAWARA-C:
        calibrationTool = import_MIAWARAC_calibrationTool(calibrationTool);
     elseif strcmp(instrumentName,'WIRAC')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        
        % Filtering options of level 1a
        calibrationTool.filterByTransmittance = false;
        calibrationTool.transmittanceThreshold = 0.05;
        calibrationTool.filterByFlags = true;
        
        %%% Flags level 1b:
        % Minimum number of averaged spectra needed for the level 1b
        calibrationTool.minNumberOfAvgSpectra = 3;
        % transmittance threshold for flagging level 1b:
        calibrationTool.troposphericTransmittanceFlag = 0.15;

        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % Import specific parameter and functions for this instrument
        calibrationTool = import_WIRAC_calibrationTool(calibrationTool);        
        % an extra folder where we copy missing anetz data (from STARTWAVE)
        calibrationTool.meteoAnetzExtraFolder = '/storage/atmosphere/MeteoSchweiz/extra/';

        % Specific to WIRAC:
        calibrationTool.outlierDectectionType = 'noFFT';

        calibrationTool.saveLevel1a = true;
        calibrationTool.integratedSpectraPlot = true;
        calibrationTool.saveLevel1b = true;

        calibrationTool.savePlanckIntensity = true;
        calibrationTool.check_deltaTC = true;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launching the calibration and integration processes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now, we keep the number of spectrometer as condition for separating 
% between GROSOM and mopi5.
    if calibrationTool.numberOfSpectrometer==1
        % Performing calibration
        if calibrate
            try
                calibrationTool = run_calibration(calibrationTool);
            catch ME
                warning(['Problem with the calibration of day: ', calibrationTool.dateStr]);
                disp(ME.message)
            end
        end
        % Performing integration
        if integrate
            try
                [calibrationTool, level1] = run_integration(calibrationTool);
            catch ME
                warning(['Problem with the integration of day: ', calibrationTool.dateStr]);
                disp(ME.message)
            end
        end
    else
        for m=1:length(modelFFTS)
            calibrationTool = calibrationTool.import_spectrometer(calibrationTool, modelFFTS(m));
            if calibrate
                try  
                    run_calibration(calibrationTool);
                catch ME
                    warning(['Problem with the calibration of ' calibrationTool.spectrometer ':']);
                    disp(ME.message)
                end
            end
            if integrate
                try
                    [calibrationTool, level1] = run_integration(calibrationTool);
                catch ME
                    warning('Problem with the integration:');
                    disp(ME.message)
                end
                %clearvars calibrationTool
            end 
        end
    end
    %clearvars level1
    clear functions   
end
