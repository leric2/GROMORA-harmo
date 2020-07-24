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
instrumentName='SOMORA';

% Type of calibration to do: standard or debug
calibrationType='standard';

calibrate = false;
integrate = true;
readLabviewLog = false;

% Define the dates for the calibration:
dates=datenum('2019_02_21','yyyy_mm_dd'):datenum('2019_02_21','yyyy_mm_dd');
% dates=[datenum('2011_03_01','yyyy_mm_dd'):datenum('2011_03_20','yyyy_mm_dd'),...
%     datenum('2011_05_01','yyyy_mm_dd'):datenum('2011_05_20','yyyy_mm_dd'),...
%     datenum('2012_03_01','yyyy_mm_dd'):datenum('2012_03_20','yyyy_mm_dd'),...
%     datenum('2012_05_01','yyyy_mm_dd'):datenum('2012_05_20','yyyy_mm_dd'),...
%     datenum('2013_03_01','yyyy_mm_dd'):datenum('2013_03_20','yyyy_mm_dd'),...
%     datenum('2013_05_01','yyyy_mm_dd'):datenum('2013_05_20','yyyy_mm_dd'),...
%     datenum('2014_03_01','yyyy_mm_dd'):datenum('2014_03_20','yyyy_mm_dd'),...
%     datenum('2014_05_01','yyyy_mm_dd'):datenum('2014_05_20','yyyy_mm_dd'),...
%     datenum('2015_03_01','yyyy_mm_dd'):datenum('2015_03_20','yyyy_mm_dd'),...
%     datenum('2015_05_01','yyyy_mm_dd'):datenum('2015_05_20','yyyy_mm_dd'),...
%     datenum('2016_03_01','yyyy_mm_dd'):datenum('2016_03_20','yyyy_mm_dd'),...
%     datenum('2016_05_01','yyyy_mm_dd'):datenum('2016_05_20','yyyy_mm_dd'),...
%     datenum('2017_03_01','yyyy_mm_dd'):datenum('2017_03_20','yyyy_mm_dd'),...
%     datenum('2017_05_01','yyyy_mm_dd'):datenum('2017_05_20','yyyy_mm_dd'),...
%     datenum('2018_01_01','yyyy_mm_dd'):datenum('2018_04_20','yyyy_mm_dd')];



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
        
        calibrationTool.filterByTransmittance = true;
    
        % Temperature of the cold load
        calibrationTool.TCold=80;
    elseif strcmp(instrumentName, 'SOMORA')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
    
        calibrationTool.filterByTransmittance = true;
        
        % Temperature of the cold load
        calibrationTool.TCold=80;        

    elseif strcmp(instrumentName,'mopi5')
        % FOR MOPI:
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
        
        calibrationTool.filterByTransmittance = true;
    
        % Temperature of the cold load
        calibrationTool.TCold=80;
        
        % the number of the spectrometer models we are interested in
        % see order in calibrationTool.spectrometerTypes
        %modelFFTS=[1 3 4];
        modelFFTS=[1 3 4 ];
        
    elseif strcmp(instrumentName,'MIAWARA-C')
        % FOR MIAWARA-C:
        % Everything stored into "import_default_calibrationTool"
    end
    
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
                calibrationTool = run_integration(calibrationTool);
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
            try
                calibrationTool = run_integration(calibrationTool);
            catch ME
                warning('Problem with the integration:');
                disp(ME.message)
            end 
        end
        
        % Integrate the successful calibration into 1 output file... maybe
        % not yet

    end
end

