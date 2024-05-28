function calibrationTool = import_MOPI5_calibrationTool(calibrationTool)
%==========================================================================
% NAME      | import_MOPI5_calibrationTool.m
% TYPE      | Function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 12.2020
%           |
% ABSTRACT  | Complete the calibrationTool structure for MOPI5
%           |
% ARGUMENTS | INPUTS: 	1. calibrationTool: the default toolbox
%           |
%           | OUTPUTS: 	2. calibrationTool: the default toolbox completer for MOPI5
%           |
% COMMENTS  | External documentation for this structure is available on the
%           | git server of IAP
%           |
%==========================================================================

% Check that instrument name corresponds to this function
assert(strcmp(calibrationTool.instrumentName, 'mopi5'),'Wrong instrument toolbox !')

calibrationTool.binaryDataExtension = '.bin';
calibrationTool.logFileDataExtension = '.txt';

calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';

calibrationTool.positionIndAsName = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.dataLocation='BERN';
calibrationTool.PI_NAME='Murk;Axel';
calibrationTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
calibrationTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
calibrationTool.PI_EMAIL='axel.murk@iap.unibe.ch';
calibrationTool.dataSource='';

% Geolocation
calibrationTool.lon=7.44;
calibrationTool.lat=46.95;
calibrationTool.altitude=560;
calibrationTool.azimuthAngle=45;

calibrationTool.timeZone = 'Z';
calibrationTool.dateTime.TimeZone = calibrationTool.timeZone;

% Observation frequency
calibrationTool.observationFreq=1.10836e11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrometer data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See the import spectrometer functions
calibrationTool.import_spectrometer = @(calibrationTool, modelFFTS) import_spectrometer_mopi5(calibrationTool, modelFFTS);
calibrationTool.numberOfSpectrometer=4;

calibrationTool.numberOfChannels=16384;

% This one should correspond to the DC channel
%calibrationTool.LOFreqTot=1.10836e11;
calibrationTool.LOFreq1 = 1.06986e11; %TODO
calibrationTool.LOFreq2 = 3.3e9;
calibrationTool.LOFreq3 = 3.8e9;
calibrationTool.DCChannel=1; %=Nchannel/2 ??

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folder, Raw and log file data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.rawFileFolder=['/storage/lake/instrumentdata/mopi5/' calibrationTool.dateStr(1:4) '/'];
%calibrationTool.rawFileFolder=['smb://resstore.unibe.ch/iap_mwlake/instrumentdata/mopi5/' calibrationTool.dateStr(1:4) '/'];

calibrationTool.level1Folder='/storage/atmosphere/instruments/mopi5/level1/';
calibrationTool.extraFileFolder='/storage/atmosphere/instruments/mopi5/extraRawFiles/'; % no write permission on the IAP lake
%calibrationTool.level1Folder='/home/eric/Documents/PhD/DATA/MOPI5/';

calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.instrumentName,'_', calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];

calibrationTool.filename=[calibrationTool.instrumentName,'_', calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];


% Log
calibrationTool.delimiter_logfile = ';';
calibrationTool.harmonize_log=@(calibrationTool, log) harmonize_log_mopi5(calibrationTool, log);
calibrationTool.THotUnit='K';

calibrationTool.indiceCold=2;
calibrationTool.indiceAntenna=5;
calibrationTool.indiceHot=1;

calibrationTool.elevationAngleAntenna=140;
calibrationTool.elevationAngleCold=265;
calibrationTool.elevationAngleHot=85;

calibrationTool.cycleDurationCold = 2;
calibrationTool.cycleDurationSky = 2;
calibrationTool.cycleDurationHot = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw File and Log

%calibrationTool.numberOfCyclesExpected=3940;
%calibrationTool.toleranceNumberCycles=0.01*calibrationTool.numberOfCyclesExpected;
calibrationTool.checkLevel0=true;
calibrationTool.tippingSize=27;
calibrationTool.flipped_spectra=false;

calibrationTool.elevationAngleTolerance=2;
calibrationTool.elevationAngleHotTol = 0;
calibrationTool.elevationAngleColdTol = 0;

%Temperature
calibrationTool.THotTh=292.2;
calibrationTool.THotAbsThresh=5;
calibrationTool.hotTemperatureStdThreshold=0.2;

calibrationTool.numberOfTippingCurveExpected=4;
calibrationTool.toleranceTippingCurves=2;

% Calibration
% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid
calibrationTool.elevationAngleTolerance=10;
calibrationTool.adcOverloadThresh = 0;
calibrationTool.stdAntAngleThresh = 0.5;

calibrationTool.minNumberOfIndicePerCycle=15;

calibrationTool.threshNumRawSpectraHot=0.1*calibrationTool.numberOfChannels;
calibrationTool.threshNumRawSpectraCold=0.1*calibrationTool.numberOfChannels;
calibrationTool.threshNumRawSpectraAnt = 0.1*calibrationTool.numberOfChannels;

calibrationTool.maxProportionOfIndLN2LevelOutlier = 0.3;
calibrationTool.maxProportionOfIndLN2SensorOutlier = 0.3;

% Filters for flagging "bad channels"
% On 10 minutes spectra
calibrationTool.maxStdDevTbCal = nan; %TODO
calibrationTool.maxStdDevTbInt = 15;

calibrationTool.filterTypeChannelQualityCal = 3;
calibrationTool.filterTypeChannelQualityInt = 4;

calibrationTool.filter1.TbMax=300;
calibrationTool.filter1.TbMin=20;
calibrationTool.filter1.boxCarSize=51;
calibrationTool.filter1.boxCarThresh=7;

% On hourly spectra
calibrationTool.filter2.TbMax=300;
calibrationTool.filter2.TbMin=20;
calibrationTool.filter2.boxCarSize=51;
calibrationTool.filter2.boxCarThresh=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read meteo data
calibrationTool.meteoFolder=['/storage/lake/instrumentdata/meteo/exwi/meteo/' calibrationTool.dateStr(1:4) '/'];

%calibrationTool.meteoFolder='smb://resstore.unibe.ch/iap_mwlake/instrumentdata/meteo/exwi/meteo/';
%calibrationTool.meteoFolder='/home/eric/Documents/PhD/METEO_DATA/';

% Function specific to this instrument
% meteo Data
%calibrationTool.get_meteo_data = @(calibrationTool,correctedSpectra) get_meteo_data_payerne(calibrationTool,correctedSpectra);
calibrationTool.read_meteo_data =@(calibrationTool) read_meteo_data_unibe(calibrationTool);
calibrationTool.add_meteo_data = @(calibrationTool, meteoData, correctedSpectra) add_meteo_data_generic(calibrationTool, meteoData, correctedSpectra);

% Backup reading of MCH ground station data (ANETZ)
calibrationTool.meteoAnetzFolder = ['/storage/atmosphere/MeteoSchweiz/' calibrationTool.dateStr(1:4) '/'];
calibrationTool.anetzStnName = 'BER';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tipping curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO
calibrationTool.doTippingCurve = false;
%calibrationTool.run_tipping_curve = @(rawSpectra, log, calibrationTool) run_tipping_curve_generic(rawSpectra,log, calibrationTool);
%calibrationTool.get_tipping_curve_data = @(rawSpectra, log, calibrationTool) get_tipping_curve_data_gromos(rawSpectra,log, calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrections
calibrationTool.transmittanceWindow=0.99;

calibrationTool.troposphericCorrection.type = 'Ingold_v1_fit';
calibrationTool.troposphericCorrection.useWings = 'both';
calibrationTool.troposphericCorrection.numberOfChannelsTropCorr = 100;
calibrationTool.troposphericCorrection.skipFraction = 0.05;
calibrationTool.troposphericCorrection.deltaT = 10.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level0 -> Level1a functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the functions that will be used for processing this
% calibration

% Selecting the functions that will be used for processing this retrieval
% Reading routine to use for the raw data
calibrationTool.read_level0=@(calibrationTool, readRawFile) read_level0_mopi5(calibrationTool, readRawFile);

% Quality check for the raw data
calibrationTool.check_level0=@(log,rawSpectra,calibrationTool) check_level0_mopi5(log,rawSpectra,calibrationTool);

% Reformatting of the raw spectra into a matrix (numberOfSpectra x
% numberOfChannels)
calibrationTool.reformat_spectra=@(rawSpectra,log,calibrationTool) reformat_spectra_generic(rawSpectra,log,calibrationTool);

% Plotting some raw spectra:
calibrationTool.plot_raw_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_raw_spectra_generic(rawSpectra,lowerLim,upperLim,N);

% TODO
% Find the sky temperature at zenith with a tipping curve
%    calibrationTool.find_T_sky_with_tipping_curve=@(rawSpectra,log,calibrationTool,calType) find_T_sky_with_tipping_curve_generic()

% Function to use for doing the calibration:
calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,calType) calibrate_generic(rawSpectra,log,calibrationTool,calType);
%calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,TCold,calType) calibrate_mopi5(rawSpectra,log,calibrationTool,TCold,calType);

% Plot some calibrated spectra:
%calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N) plot_spectra_mopi(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N);
calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,N) plot_spectra_generic(calibrationTool,drift,meteoData, calibratedSpectra,N);

% Function for quality check of the calibrated spectra
%calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_generic(log,calibrationTool,calibratedSpectra);
calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_mopi5(log,calibrationTool,calibratedSpectra);

% Function saving the calibrated spectra into netCDF file
calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(calibrationTool,log,calibratedSpectra,warningLevel0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level1a -> Level1b functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function reading the daily calibrated spectra from netCDF file
calibrationTool.read_level1a = @(calibrationTool) read_level1a_daily(calibrationTool);

% Check of the channels quality on the calibrated spectra:
calibrationTool.check_channel_quality= @(calibratedSpectra,calibrationTool,filterN) check_channel_quality_generic(calibratedSpectra,calibrationTool,filterN);

% Integration of level1a data
calibrationTool.integrate_calibrated_spectra= @(calibrationTool,calibratedSpectra) integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra);

% Function for plotting the integrated spectra (when hourly)
calibrationTool.plot_integrated_spectra = @(calibrationTool,rawSpectra) plot_integrated_spectra_generic(calibrationTool,rawSpectra);

calibrationTool.tropospheric_correction = @(integration,calibrationTool) tropospheric_correction_generic(integration,calibrationTool);

% Window correction for the calibrated spectra
calibrationTool.window_correction= @(calibrationTool,level1b) window_correction_generic(calibrationTool,level1b);

% Check of the integrated spectra
calibrationTool.check_integrated = @(calibrationTool,integratedSpectra) check_integrated_generic(calibrationTool,integratedSpectra);

% Function saving the calibrated spectra into netCDF file
calibrationTool.save_level1b=@(calibrationTool,level1b) save_level1b_daily(calibrationTool,level1b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter varying with time for the instruments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calibrationTool.timeNumber<datenum(2019,02,12)
    calibrationTool.elevationAngleHot=90;
end