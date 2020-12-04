function calibrationTool = import_MIAWARAC_calibrationTool(calibrationTool)
%==========================================================================
% NAME      | import_MIAWARAC_calibrationTool(calibrationTool)
% TYPE      | Function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 12.2020
%           |
% ARGUMENTS | INPUTS: - calibrationTool: the default toolbox
%           |
%           | OUTPUTS: - calibrationTool: the default toolbox for MIAWARA-C
%           |
%==========================================================================

% Check that instrument name corresponds to this function
assert(strcmp(calibrationTool.instrumentName, 'MIAWARA-C'),'Wrong instrument toolbox !')

calibrationTool.binaryDataExtension = '.bin';
calibrationTool.logFileDataExtension = '.txt';

calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';

calibrationTool.positionIndAsName = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta data
calibrationTool.dataLocation='NY.ALESUND';
calibrationTool.PI_NAME='Murk;Axel';
calibrationTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
calibrationTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
calibrationTool.PI_EMAIL='axel.murk@iap.unibe.ch';
calibrationTool.dataSource='MWR.H2O_UBERN002';

% Geolocation
calibrationTool.lon=11.92;
calibrationTool.lat=78.92;
calibrationTool.altitude=0;
calibrationTool.azimuthAngle=0;


% Observation frequency
%         calibrationTool.observationFreq=%;
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Spectrometer data:
calibrationTool.numberOfSpectrometer=1;
calibrationTool.spectrometer='AC240';
%
%         calibrationTool.samplingRateFFTS=%2000;
%
calibrationTool.numberOfChannels=32768;
%
%         % Local oscillators information
%         calibrationTool.fLO1=%1.45875e11;
%         calibrationTool.fLO2=%3.6e9;
%
%         % This one should correspond to the DC channel
%         calibrationTool.LOFreqTot=%calibrationTool.fLO1-calibrationTool.fLO2-0.5e9;
%
%         %calibrationTool.DCChannel=16384;
%         calibrationTool.instrumentBandwidth=%1e9;
%         %calibrationTool.LOFreq=1.45875e11;
%
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Parameters for flagging
%         calibrationTool.TSysCenterTh=%2750;
%         calibrationTool.TSysThresh=%100;
%         calibrationTool.stdTSysThresh=%8;
%
%         calibrationTool.THotTh=%313.9;
%         calibrationTool.THotAbsThresh=%2;
%
%
%         calibrationTool.numberOfTippingCurveExpected=%48;
%         calibrationTool.toleranceTippingCurves=%2;
%         calibrationTool.elevationAngleAntenna=%40;
%         calibrationTool.elevationAngleCold=%-84;
%         % -84 before 2020 ?????
%         %calibrationTool.elevationAngleCold=-84;
%         calibrationTool.elevationAngleHot=%160;
%
%         calibrationTool.elevationAngleTolerance=%5;
%         % Considering the expected number of tipping curve:
%         calibrationTool.numberOfCyclesExpected=%1500;
%         calibrationTool.toleranceNumberCycles=%0.01*calibrationTool.numberOfCyclesExpected;
%         calibrationTool.tippingSize=%27;
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Corrections
%         calibrationTool.tWindow=%0.99;
%
calibrationTool.checkLevel0=false;
%
calibrationTool.flipped_spectra=false;
%         calibrationTool.flip_spectra=%@(rawSpectra) flip_spectra_gromos(rawSpectra);
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Log file from the instrument:
%         calibrationTool.THotUnit='degreeC';

%         % Function for the harmonization of the log
calibrationTool.harmonize_log=@(calibrationTool, log) harmonize_log_miawarac(calibrationTool, log);
%
%         calibrationTool.indiceCold=%0;
%         calibrationTool.indiceAntenna=%1;
%         calibrationTool.indiceHot=%2;
%
%         % Calibration outlier management
%         calibrationTool.threshNumRawSpectraHot=%0.05*calibrationTool.numberOfChannels;
%         calibrationTool.threshNumRawSpectraCold=%0.05*calibrationTool.numberOfChannels;
calibrationTool.delimiter_logfile = ';';
% paths
%calibrationTool.rawFileFolder='/home/franziska/Documents/MW/play_MIA-C_calibration/';%['/mnt/instrumentdata/miawarac/' dateStr(1:4) '/'];
calibrationTool.extraFileFolder='/scratch/GROSOM/ExtraRawFiles/'; % no write permission on the IAP lake
calibrationTool.rawFileFolder=['/scratch/'];
calibrationTool.level1Folder='/home/franziska/Documents/MW/play_MIA-C_calibration/';

calibrationTool.filename=[calibrationTool.instrumentName,'_', calibrationTool.dateStr(1:4) '_' calibrationTool.dateStr(6:7) '_' calibrationTool.dateStr(9:10)];
calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];

calibrationTool.meteoFolder='/home/franziska/Documents/MW/play_MIA-C_calibration/';
calibrationTool.observationFreq=22.235;

calibrationTool.calibrationTime=60;

%calibrationTool.fLO1=1.49275e11;
%calibrationTool.fLO2=5.6e9;
%calibrationTool.fLO3=2e9;

% This one should correspond to the DC channel
calibrationTool.LOFreqTot=1.10e11;
calibrationTool.DCChannel=1; %=Nchannel/2 ??

% files for the calibration
calibrationTool.channel_freqs = 'miawarac_mysql_channels_freqs.dat';
calibrationTool.elcorr_file   = 'miawarac_elcorr.mat';
calibrationTool.antenna_file  = 'miawarac_antenna.txt';

calibrationTool.read_level0=@(calibrationTool, rawFileReading) read_level0_missing(calibrationTool, rawFileReading);

calibrationTool.filenameLevel1a=['/home/franziska/Documents/MW/play_MIA-C_calibration/MIAWARA-C_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];

calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,calType) run_balancing_calibration(rawSpectra,log,calibrationTool,calType);
calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_miawara_c(log,calibrationTool,calibratedSpectra);

% tipping curve
calibrationTool.run_tipping_curve = @(rawSpectra, log, calibrationTool) run_tipping_curve_generic(rawSpectra,log, calibrationTool);
calibrationTool.get_tipping_curve_data = @(rawSpectra, log, calibrationTool) get_tipping_curve_data_miawarac(rawSpectra,log, calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From calibration GROSOM
% Type of calibration to do: standard of debug
calibrationType='standard';

% working directory
root_dir = '/home/franziska/Documents/MW/GROSOM-harmo/';
%cd work_path
addpath(genpath(root_dir))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the calibrationTool for this particular day and instrument:
% calibrationTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};

% for gaining time.
calibrationTool.level1aExist=false;

% Time interval for doing the calibration
calibrationTool.calibrationTime=360;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter varying with time for the instruments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Meteo
YYYY = calibrationTool.dateStr(1:4);
MM   = calibrationTool.dateStr(6:7);
DD   = calibrationTool.dateStr(9:10);

if calibrationTool.timeNumber > datenum('2015-09-01') && calibrationTool.timeNumber < datenum(2016,11,25)
    calibrationTool.get_meteo_data = @(calibrationTool) get_meteo_data_grc(calibrationTool);
    %calibrationTool.meteoFile = ['/data/miradara3/instrumentdata/gromosc/' YYYY '/GROMOS-C*' YYYY '_' MM '_' DD];
    calibrationTool.meteoFile = [calibrationTool.rawFileFolder 'GROMOS-C*' YYYY '_' MM '_' DD];
elseif calibrationTool.timeNumber > datenum('2017-04-06 23:00')
    calibrationTool.get_meteo_data = @(calibrationTool) get_meteo_data_miac_ownfile(calibrationTool);
    calibrationTool.meteoFile = [ calibrationTool.rawFileFolder 'MIAWARA-C_meteo_' YYYY MM DD];
else
    calibrationTool.get_meteo_data = @(log) get_meteo_data_miac(log);
end
