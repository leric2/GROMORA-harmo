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
%           | OUTPUTS: Level1a and level1b, saved from the run_calibration 
%           | function
%           |
% CALLS     | import_default_calibrationTool(instrumentName,dateStr)
%           | run_calibration(calibrationTool)
%           |
%           |
%==========================================================================htop


clear all; close all; clc; clear functions; %clear mex;

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C' // 'WIRAC' //
instrumentName='GROMOS';

% Type of calibration to do: standard or debug
calibrationType='standard';

calibrate = false;
integrate = true;
readLabviewLog = false;

% GROMOS from 10.03.2010 only (after change in SW, see logfile), meteo from
% 12.05.2010
warning off backtrace
% Define the dates for the calibration:
%dates=datenum('2009_01_01','yyyy_mm_dd'):datenum('2009_01_01','yyyy_mm_dd');
dates=datenum('1995_05_15','yyyy_mm_dd'):datenum('1995_05_15','yyyy_mm_dd');
% dates=[datenum('2009_08_22','yyyy_mm_dd'):datenum('2009_08_23','yyyy_mm_dd'),...
%      datenum('2009_09_21','yyyy_mm_dd'):datenum('2009_09_21','yyyy_mm_dd')];
%      datenum('2010_08_02','yyyy_mm_dd'):datenum('2010_08_02','yyyy_mm_dd'),...
%      datenum('2010_08_31','yyyy_mm_dd'):datenum('2010_08_31','yyyy_mm_dd')];
%      datenum('2020_10_10','yyyy_mm_dd'):datenum('2020_10_19','yyyy_mm_dd'),...
%       datenum('2020_11_18','yyyy_mm_dd'):datenum('2020_11_19','yyyy_mm_dd'),...
%       datenum('2020_12_18','yyyy_mm_dd'):datenum('2020_12_19','yyyy_mm_dd'),...
%       datenum('2020_01_18','yyyy_mm_dd'):datenum('2020_01_19','yyyy_mm_dd'),...
%       datenum('2020_02_18','yyyy_mm_dd'):datenum('2020_02_19','yyyy_mm_dd'),...
%       datenum('2020_03_18','yyyy_mm_dd'):datenum('2020_03_19','yyyy_mm_dd'),...
%       datenum('2020_04_18','yyyy_mm_dd'):datenum('2020_04_19','yyyy_mm_dd'),...
%       datenum('2020_05_18','yyyy_mm_dd'):datenum('2020_05_19','yyyy_mm_dd'),...
%       datenum('2020_06_18','yyyy_mm_dd'):datenum('2020_06_19','yyyy_mm_dd')];

% good_date mopi5
% dates=[datenum('2019_01_03','yyyy_mm_dd'):datenum('2019_01_09','yyyy_mm_dd'),...
%     datenum('2019_01_30','yyyy_mm_dd'):datenum('2019_02_22','yyyy_mm_dd'),...
%     datenum('2019_03_01','yyyy_mm_dd'):datenum('2019_03_01','yyyy_mm_dd'),...
%     datenum('2019_03_12','yyyy_mm_dd'):datenum('2019_03_12','yyyy_mm_dd'),...
%     datenum('2019_04_25','yyyy_mm_dd'):datenum('2019_05_04','yyyy_mm_dd'),...
%     datenum('2019_06_11','yyyy_mm_dd'):datenum('2019_06_18','yyyy_mm_dd')];


if (strcmp(instrumentName,'GROMOS') | strcmp(instrumentName,'SOMORA')) & readLabviewLog
    labviewLogFolder = '/home/esauvageat/Documents/GROMORA/Analysis/InputsCalibration/';
    labviewLogFolder = '/storage/tub/instruments/gromos/level1/GROMORA/InputsCalibration/';
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
    
    calibrationTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
    calibrationTool.referenceTime = datenum(1970,1,1,0,0,0);
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
    calibrationTool.saveLevel1a = false;
    calibrationTool.integratedSpectraPlot = true;
    calibrationTool.saveLevel1b = true;
    
    calibrationTool.savePlanckIntensity = false;
    calibrationTool.check_deltaTC = false;
    
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
        calibrationTool = import_GROMOS_FB_calibrationTool(calibrationTool);

        % an extra folder where we copy missing anetz data (from STARTWAVE)
        calibrationTool.meteoAnetzExtraFolder = '/storage/tub/MeteoSchweiz/extra/';

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
        calibrationTool.meteoAnetzExtraFolder = '/storage/tub/MeteoSchweiz/extra/';

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
        calibrationTool.meteoAnetzExtraFolder = '/storage/tub/MeteoSchweiz/extra/';

        % Specific to WIRAC:
        calibrationTool.outlierDectectionType = 'noFFT';

        calibrationTool.saveLevel1a = true;
        calibrationTool.integratedSpectraPlot = true;
        calibrationTool.saveLevel1b = true;

        calibrationTool.savePlanckIntensity = true;
        calibrationTool.check_deltaTC = true;
    end
    
    %calibrationTool.transmittanceWindow  = calibrationTool.transmittanceWindow *0.97
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launching the calibration and integration processes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now, we keep the number of spectrometer as condition for separating 
% between GROSOM and mopi5.
        % Performing integration
        if integrate
            try
                [calibrationTool, level1] = run_integration(calibrationTool);
            catch ME
                warning('Problem integrating this day !')
                disp(ME.message)
            end  
        end

    %clearvars level1
    clear functions   
end

