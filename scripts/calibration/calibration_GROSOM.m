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

clear; close all; clc; clear functions; %clear mex;

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C'
instrumentName='GROMOS';

% Type of calibration to do: standard or debug
calibrationType='standard';

calibrate = false;
integrate = false;
readLabviewLog = false;

% GROMOS from 10.03.2010 only (after change in SW, see logfile), meteo from
% 12.05.2010

% Define the dates for the calibration:
dates=datenum('2019_01_04','yyyy_mm_dd'):datenum('2019_01_04','yyyy_mm_dd');

% good_date mopi5
% dates=[datenum('2019_01_03','yyyy_mm_dd'):datenum('2019_01_09','yyyy_mm_dd'),...
%     datenum('2019_01_30','yyyy_mm_dd'):datenum('2019_02_22','yyyy_mm_dd'),...
%     datenum('2019_03_01','yyyy_mm_dd'):datenum('2019_03_01','yyyy_mm_dd'),...
%     datenum('2019_03_12','yyyy_mm_dd'):datenum('2019_03_12','yyyy_mm_dd'),...
%     datenum('2019_04_25','yyyy_mm_dd'):datenum('2019_05_04','yyyy_mm_dd'),...
%     datenum('2019_06_11','yyyy_mm_dd'):datenum('2019_06_18','yyyy_mm_dd')];

% dates=[datenum('2019_11_01','yyyy_mm_dd'):datenum('2019_11_01','yyyy_mm_dd'),...
%     datenum('2009_10_01','yyyy_mm_dd'):datenum('2009_10_21','yyyy_mm_dd'),...
%     datenum('2010_04_01','yyyy_mm_dd'):datenum('2010_04_21','yyyy_mm_dd'),...
%     datenum('2010_10_01','yyyy_mm_dd'):datenum('2010_10_21','yyyy_mm_dd'),...
%     datenum('2011_04_01','yyyy_mm_dd'):datenum('2011_04_21','yyyy_mm_dd'),...
%     datenum('2011_10_01','yyyy_mm_dd'):datenum('2011_10_21','yyyy_mm_dd'),...
%     datenum('2012_04_01','yyyy_mm_dd'):datenum('2012_04_21','yyyy_mm_dd'),...
%     datenum('2012_10_01','yyyy_mm_dd'):datenum('2012_10_21','yyyy_mm_dd'),...
%     datenum('2013_04_01','yyyy_mm_dd'):datenum('2013_04_21','yyyy_mm_dd'),...
%     datenum('2013_10_01','yyyy_mm_dd'):datenum('2013_10_21','yyyy_mm_dd'),...
%     datenum('2014_04_01','yyyy_mm_dd'):datenum('2014_04_21','yyyy_mm_dd'),...
%     datenum('2014_10_01','yyyy_mm_dd'):datenum('2014_10_21','yyyy_mm_dd'),...
%     datenum('2015_04_01','yyyy_mm_dd'):datenum('2015_04_21','yyyy_mm_dd'),...
%     datenum('2015_10_01','yyyy_mm_dd'):datenum('2015_10_21','yyyy_mm_dd'),...
%     datenum('2016_04_01','yyyy_mm_dd'):datenum('2016_04_21','yyyy_mm_dd'),...
%     datenum('2016_10_01','yyyy_mm_dd'):datenum('2016_10_21','yyyy_mm_dd'),...
%     datenum('2017_04_01','yyyy_mm_dd'):datenum('2017_04_21','yyyy_mm_dd'),...
%     datenum('2017_10_01','yyyy_mm_dd'):datenum('2017_10_21','yyyy_mm_dd'),...
%     datenum('2018_04_01','yyyy_mm_dd'):datenum('2018_04_21','yyyy_mm_dd'),...
%     datenum('2018_10_01','yyyy_mm_dd'):datenum('2018_10_21','yyyy_mm_dd'),...
%     datenum('2019_04_01','yyyy_mm_dd'):datenum('2019_04_21','yyyy_mm_dd'),...
%     datenum('2019_10_01','yyyy_mm_dd'):datenum('2019_10_21','yyyy_mm_dd'),...
%     datenum('2020_04_01','yyyy_mm_dd'):datenum('2020_04_21','yyyy_mm_dd')];

%dates=datenum('2015_09_27','yyyy_mm_dd')
if (strcmp(instrumentName,'GROMOS') | strcmp(instrumentName,'SOMORA')) & readLabviewLog
    labviewLogFolder = '/home/esauvageat/Documents/GROSOM/Analysis/InputsCalibration/';
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
    calibrationTool.calibrationVersion = '1.0';
    
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
    
    calibrationTool.rawSpectraPlot = false;
    calibrationTool.calibratedSpectraPlot = true;
    calibrationTool.calibratedSpectraSpectralPlot = true;
    calibrationTool.saveLevel1a = true;
    calibrationTool.integratedSpectraPlot = true;
    calibrationTool.saveLevel1b = true;
    
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
        calibrationTool.transmittanceThreshold = 0.2;
        calibrationTool.filterByFlags = true;
        
        % Minimum number of averaged spectra needed for the level 1b
        calibrationTool.minNumberOfAvgSpectra = 3;

        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % Import specific parameter and functions for this instrument
        calibrationTool = import_GROMOS_calibrationTool(calibrationTool);
        
    elseif strcmp(instrumentName, 'SOMORA')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        
        % Filtering options of level 1a
        calibrationTool.filterByTransmittance = true;
        calibrationTool.transmittanceThreshold = 0.2;
        calibrationTool.filterByFlags = true;
        
        % Minimum number of averaged spectra needed for the level 1b
        calibrationTool.minNumberOfAvgSpectra = 3;

        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % Import specific parameter and functions for this instrument
        calibrationTool = import_SOMORA_calibrationTool(calibrationTool);

    elseif strcmp(instrumentName,'mopi5')
        % FOR MOPI:
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=6*60;
        calibrationTool.minNumberOfAvgSpectra = 6;
        
        calibrationTool.filterByTransmittance = false; % Best to keep false for MOPI5 studies
        calibrationTool.transmittanceThreshold = 0.2;
        calibrationTool.filterByFlags = false;
        
        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % the number of the spectrometer models we are interested in
        % see order in calibrationTool.spectrometerTypes
        modelFFTS=[1 3 4];
        
        calibrationTool = import_MOPI5_calibrationTool(calibrationTool);
    elseif strcmp(instrumentName,'MIAWARA-C')
        % FOR MIAWARA-C:
        calibrationTool = import_MIAWARAC_calibrationTool(calibrationTool);
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
                [calibrationTool, level1b] = run_integration(calibrationTool);
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
                    [calibrationTool, level1b] = run_integration(calibrationTool);
                catch ME
                    warning('Problem with the integration:');
                    disp(ME.message)
                end
            end 
        end
    end
    % clearvars level1b, calibrationTool
    clear functions   
end

