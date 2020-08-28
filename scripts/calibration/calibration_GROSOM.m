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

% Type of calibration to do: standard or debug
calibrationType='standard';

calibrate = true;
integrate = true;
readLabviewLog = true;

% GROMOS from 10.03.2010 only (after change in SW, see logfile), meteo from
% 12.05.2010

% Define the dates for the calibration:
dates=datenum('2019_01_01','yyyy_mm_dd'):datenum('2019_01_02','yyyy_mm_dd');
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
    labviewLog = read_labview_log_generic(instrumentName);
else
    labviewLog = struct();
end

% working directory
%root_dir = '/home/franziska/Documents/MW/GROSOM-harmo/';
root_dir = pwd;
%cd work_path
addpath(genpath(root_dir))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining all parameters for the calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for d = 1:numel(dates)
    dateStr=datestr(dates(d),'yyyy_mm_dd');
    
    % Import default tools for running a retrieval for a given instrument
    calibrationTool=import_default_calibrationTool(instrumentName,dateStr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing the calibrationTool for this particular day and instrument:
    % calibrationTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};
    
    % for gaining time.
    %calibrationTool.level1aExist=false;
    
    calibrationTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
    calibrationTool.calendar='standard';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    calibrationTool.hotSpectraNumberOfStdDev=3;
    calibrationTool.coldSpectraNumberOfStdDev=3;
    
    calibrationTool.labviewLog = labviewLog;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug mode and plot options
    calibrationTool.calType=calibrationType;
    
    calibrationTool.rawSpectraPlot=false;
    calibrationTool.calibratedSpectraPlot=true;
    calibrationTool.integratedSpectraPlot=true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instrument specific parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % On the long term this should be all taken from import_default_calTool
    
    if strcmp(instrumentName,'GROMOS')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        calibrationTool.minNumberOfAvgSpectra = 3;
        
        calibrationTool.filterByTransmittance = true;
        calibrationTool.filterByFlags = true;
    
        % Temperature of the cold load
        calibrationTool.TCold=80;
    elseif strcmp(instrumentName, 'SOMORA')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        calibrationTool.minNumberOfAvgSpectra = 3;
    
        calibrationTool.filterByTransmittance = true;
        calibrationTool.filterByFlags = true;
        
        % Temperature of the cold load
        calibrationTool.TCold=80;        

    elseif strcmp(instrumentName,'mopi5')
        % FOR MOPI:
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        calibrationTool.minNumberOfAvgSpectra = 12;
        
        calibrationTool.filterByTransmittance = false; % Best to keep false for MOPI5 studies
        calibrationTool.filterByFlags = true;
        
        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % the number of the spectrometer models we are interested in
        % see order in calibrationTool.spectrometerTypes
        %modelFFTS=[1 3 4];
        modelFFTS=[1 3 4];
        
    elseif strcmp(instrumentName,'MIAWARA-C')
        % FOR MIAWARA-C:
        % Everything stored into "import_default_calibrationTool"
    end
    
    %TODO
    % create all filenames at that stage and store them in calibrationTool
    % Not possible before if we want to include the calibrationTime for
    % instance.
    % calibrationTool = create_filenames(calibrationTool);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launching the calibration process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now, we keep it dirty for separating between GROSOM and MOPI
    if calibrationTool.numberOfSpectrometer==1
            % if commented, nothing happens --> developping purposes
        if calibrate
            try
                calibrationTool = run_calibration(calibrationTool);
            catch ME
                warning('Problem with the calibration:');
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
        
        % Integrate the successful calibration into 1 output file... maybe
        % not yet

    end
    clearvars level1b
end

