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

% 'GROMOS' // 'SOMORA' // 'MOPI5' // 'MIAWARA-C'
instrumentName='GROMOS';

% Type of calibration to do: standard of debug
calibrationType='standard';

% Define the dates for the calibration:
dates=datenum('2019_01_11','yyyy_mm_dd'):datenum('2019_01_11','yyyy_mm_dd');
%dates=datenum('2015_09_27','yyyy_mm_dd')

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
    calibrationTool.level1aExist=false;
    
    calibrationTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
    calibrationTool.calendar='standard';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    calibrationTool.hotSpectraNumberOfStdDev=3;
    calibrationTool.coldSpectraNumberOfStdDev=3;
    
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
    
        % Temperature of the cold load
        calibrationTool.TCold=80;
    elseif strcmp(instrumentName, 'SOMORA')
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
    
        % Temperature of the cold load
        calibrationTool.TCold=80;        
        % TOCHANGE
        calibrationTool.meteoFolder='/home/esauvageat/Documents/GROSOM/Analysis/MeteoFile/METEO_DATA/';
        %calibrationTool.meteoFolder='/home/eric/Documents/PhD/METEO_DATA/';
    
        % Function specific to this instrument
        % meteo Data
        calibrationTool.get_meteo_data = @(calibrationTool,correctedSpectra) get_meteo_data_payerne(calibrationTool,correctedSpectra);
    elseif strcmp(instrumentName,'MOPI5')
        % FOR MOPI:
        % Everything stored into "import_default_calibrationTool"
        calibrationTool.extraFileFolder='/home/esauvageat/Documents/GROSOM/Analysis/InputsCalibration/'; % no write permission on the IAP lake
        
        % Time interval for doing the calibration
        calibrationTool.calibrationTime=10;
    
        % Total integration time
        calibrationTool.integrationTime=60;
    
        % Temperature of the cold load
        calibrationTool.TCold=80;
        
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
            if ~calibrationTool.level1aExist
                calibrationTool = run_calibration(calibrationTool);
            end
            calibrationTool = run_integration(calibrationTool);
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

