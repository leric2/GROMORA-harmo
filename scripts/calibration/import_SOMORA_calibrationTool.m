function calibrationTool = import_SOMORA_calibrationTool(calibrationTool)
%==========================================================================
% NAME      | import_SOMORA_calibrationTool.m
% TYPE      | Function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 12.2020
%           |
% ABSTRACT  | Complete the calibrationTool structure for SOMORA
%           |
% ARGUMENTS | INPUTS: 	1. calibrationTool: the default toolbox
%           |
%           | OUTPUTS: 	2. calibrationTool: the default toolbox completer for SOMORA
%           |
% COMMENTS  | External documentation for this structure is available on the
%           | git server of IAP
%           |
%==========================================================================

% Check that instrument name corresponds to this function
assert(strcmp(calibrationTool.instrumentName, 'SOMORA'),'Wrong instrument toolbox !')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.dataLocation='PAYERNE';
calibrationTool.PI_NAME='';
calibrationTool.PI_AFFILIATION='Swiss Meteorological Institute;MCH';
calibrationTool.PI_ADDRESS='';
calibrationTool.PI_EMAIL='';
calibrationTool.dataSource='MWR.O3_MCH';

% Geolocation
calibrationTool.lon=6.95;
calibrationTool.lat=46.82;
calibrationTool.altitude=491;
calibrationTool.azimuthAngle=34;

calibrationTool.timeZone = 'Z';
calibrationTool.dateTime.TimeZone = calibrationTool.timeZone;

% Observation frequency
calibrationTool.observationFreq=1.4217504e11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrometer data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.numberOfSpectrometer=1;
calibrationTool.spectrometer='AC240';

calibrationTool.IQProcessing = false;

calibrationTool.fLO1=1.49275e11;
calibrationTool.fLO2=5.6e9;
calibrationTool.fLO3=2e9;

% This one should correspond to the DC channel
calibrationTool.LOFreqTot=calibrationTool.fLO1-calibrationTool.fLO2-calibrationTool.fLO3;
%calibrationTool.DCChannel=1; %=Nchannel/2 ??

calibrationTool.numberOfChannels=16384;

calibrationTool.instrumentBandwidth=1e9;
calibrationTool.LOFreq=1.4927504e11; %TOCHECK
calibrationTool.samplingRateFFTS=2000;

calibrationTool.badChannels=1:104;

calibrationTool.flipped_spectra=false;

calibrationTool.numberOfAquisitionSpectraHot=60000;
calibrationTool.numberOfAquisitionSpectraAntenna=120000;
calibrationTool.numberOfAquisitionSpectraCold=120000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folder, Raw and log file data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.binaryDataExtension = '.bin';
calibrationTool.logFileDataExtension = '.txt';

calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';

%calibrationTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' calibrationTool.dateStr(6:7) '/'];
calibrationTool.rawFileFolder=['/scratch/SOMORA_rawData/' calibrationTool.dateStr(1:4) '/' calibrationTool.dateStr(6:7) '/'];

%calibrationTool.rawFileFolder='/home/eric/Documents/PhD/GROSOM/rawData/';
calibrationTool.extraFileFolder='/scratch/GROSOM/ExtraRawFiles/'; % no write permission on the IAP lake
%calibrationTool.level1Folder='/home/eric/Documents/PhD/GROSOM/Level1/';
calibrationTool.level1Folder='/scratch/GROSOM/Level1/SOMORA/';
calibrationTool.filename=[calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];

% calibrationTool.rawFile = [calibrationTool.rawFileFolder, calibrationTool.filename, calibrationTool.binaryDataExtension];
% calibrationTool.logFile = [calibrationTool.rawFileFolder, calibrationTool.filename, calibrationTool.logFileDataExtension];


% Defining level1a filename to read (to be adapted for other users)
calibrationTool.filenameLevel1a = [calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];

% Log
calibrationTool.delimiter_logfile = '\t';
calibrationTool.harmonize_log=@(calibrationTool, log) harmonize_log_somora(calibrationTool, log);
calibrationTool.THotUnit='K';

calibrationTool.positionIndAsName = false;
calibrationTool.indiceCold=0;
calibrationTool.indiceAntenna=1;
calibrationTool.indiceHot=2;

calibrationTool.elevationAngleAntenna=38;
calibrationTool.elevationAngleCold=-90;
calibrationTool.elevationAngleHot=180;

calibrationTool.elevationAngleTolerance=2;
calibrationTool.elevationAngleHotTol = 0;
calibrationTool.elevationAngleColdTol = 0;

calibrationTool.cycleDurationCold = 4;
calibrationTool.cycleDurationSky = 2;
calibrationTool.cycleDurationHot = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw File and Log
calibrationTool.numberOfCyclesExpected=3940;
calibrationTool.toleranceNumberCycles=0.01*calibrationTool.numberOfCyclesExpected;

calibrationTool.tippingSize=5;

%Temperature
calibrationTool.TNoiseCenterTh=2750;
calibrationTool.TNoiseThresh=300;
calibrationTool.stdTNoiseThresh=8;

calibrationTool.THotTh=311.1;
calibrationTool.THotAbsThresh=5;
calibrationTool.hotTemperatureStdThreshold=0.05;
calibrationTool.checkLevel0=true;

calibrationTool.numberOfTippingCurveExpected=4;
calibrationTool.toleranceTippingCurves=2;

% Calibration
% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid

calibrationTool.stdAntAngleThresh = 0.5;

calibrationTool.minNumberOfIndicePerCycle=40;
calibrationTool.threshNumRawSpectraHot=0.05*calibrationTool.numberOfChannels;
calibrationTool.threshNumRawSpectraCold=0.05*calibrationTool.numberOfChannels;
calibrationTool.threshNumRawSpectraAnt = 0.05*calibrationTool.numberOfChannels;

calibrationTool.frequencyBandAroundCenterTNoise = 200e6;

calibrationTool.maxProportionOfIndLN2LevelOutlier = 0.3;
calibrationTool.maxProportionOfIndLN2SensorOutlier = 0.3;

% Filters for flagging "bad channels"
calibrationTool.maxStdDevTbCal = 25;
calibrationTool.maxStdDevTbInt = 10;

calibrationTool.filterTypeChannelQualityCal = 3;
calibrationTool.filterTypeChannelQualityInt = 2;

% On 10 minutes spectra
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
calibrationTool.meteoFolder=['/scratch/GROSOM/MeteoDataSOMORA/METEO_DATA_' calibrationTool.dateStr(1:4) '/'];
%calibrationTool.meteoFolder='/home/eric/Documents/PhD/GROSOM/METEO_DATA/';

% Function specific to this instrument
% meteo Data
%calibrationTool.get_meteo_data = @(calibrationTool,correctedSpectra) get_meteo_data_payerne(calibrationTool,correctedSpectra);
calibrationTool.read_meteo_data =@(calibrationTool) read_meteo_data_payerne(calibrationTool);
calibrationTool.add_meteo_data = @(calibrationTool, meteoData, correctedSpectra) add_meteo_data_generic(calibrationTool, meteoData, correctedSpectra);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tipping curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO
calibrationTool.doTippingCurve = false;
%calibrationTool.run_tipping_curve = @(rawSpectra, log, calibrationTool) run_tipping_curve_generic(rawSpectra,log, calibrationTool);
%calibrationTool.get_tipping_curve_data = @(rawSpectra, log, calibrationTool) get_tipping_curve_data_gromos(rawSpectra,log, calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrections
calibrationTool.troposphericCorrection.type = 'Ingold_v1';
calibrationTool.troposphericCorrection.useWings = 'both';
calibrationTool.troposphericCorrection.numberOfChannelsTropCorr = 50;
calibrationTool.troposphericCorrection.skipFraction = 0.05;
calibrationTool.troposphericCorrection.deltaT = 10.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level0 -> Level1a functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the functions that will be used for processing this
% calibration

% Selecting the functions that will be used for processing this retrieval
% Reading routine to use for the raw data
calibrationTool.read_level0=@(calibrationTool, readRawFile) read_level0_generic(calibrationTool, readRawFile);

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
calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,calType) calibrate_generic(rawSpectra,log,calibrationTool,calType);

% Plot some calibrated spectra:
calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,N) plot_spectra_generic(calibrationTool,drift,meteoData, calibratedSpectra,N);

% Function for quality check of the calibrated spectra
calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_generic(log,calibrationTool,calibratedSpectra);

% Function saving the calibrated spectra into netCDF file
calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(calibrationTool,log,calibratedSpectra,warningLevel0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level1a -> Level1b functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function reading the daily calibrated spectra from netCDF file
calibrationTool.read_level1a = @(calibrationTool) read_level1a_daily(calibrationTool);
%calibrationTool.read_level1a = @(calibrationTool,sublevel) read_level1_GROSOM(calibrationTool,sublevel);

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
if calibrationTool.timeNumber>datenum(2016,01,01) && calibrationTool.timeNumber<datenum(2017,04,01)
    calibrationTool.THotTh = 297.4;
end

% if calibrationTool.timeNumber<datenum(2019,01,01)
%     calibrationTool.rawFileFolder=['/media/esauvageat/INTENSO/RAW_DATA/' calibrationTool.dateStr(1:4) '/' calibrationTool.dateStr(6:7) '/'];
%     calibrationTool.filename=[calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
%     calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];
% end
% window transmission
if calibrationTool.timeNumber<datenum(2007,07,03) %attn anciennement  733957 cad 3 juillet 2009 mais FAUX
    calibrationTool.tWindow = 0.987; %value from 2002 to 3/7/2007 (att. window has been changed 8/2/2006 but new t not measured) nb: ema 9/11/2007
elseif (calibrationTool.timeNumber>= datenum(2007,07,03) && calibrationTool.timeNumber<datenum(2010,6,26)) %ie 734315
    calibrationTool.tWindow = 0.979; %value from 3/7/2007 to 25/6/2010   nb: ema 9/11/2007; calc sur la base de R-J
elseif (calibrationTool.timeNumber>= datenum(2010,06,25) && calibrationTool.timeNumber<datenum(2012,03,02)) %ie 734931 2/3/2012
    calibrationTool.tWindow = 0.985; %value from 26/6/2010 to 2/3/2012 %attn anciennement 0.976 ???
elseif (calibrationTool.timeNumber>= datenum(2012,03,02) && calibrationTool.timeNumber<datenum(2014,06,16)) %ie 16/06/2014
    calibrationTool.tWindow = 0.999; %value from 3/3/2012 to 15/06/2014
elseif (calibrationTool.timeNumber>= datenum(2014,06,16) && calibrationTool.timeNumber<datenum(2015,09,08))
    calibrationTool.tWindow = 1.0105; %value from 16/6/2014 to 07/09/2015
elseif (calibrationTool.timeNumber>= datenum(2015,09,09) && calibrationTool.timeNumber<datenum(2017,05,30))
    calibrationTool.tWindow = 1.009; %value from 8/09/2015 to 29/05/2017
elseif (calibrationTool.timeNumber>= datenum(2017,05,31) && calibrationTool.timeNumber<datenum(2018,08,21))
    calibrationTool.tWindow = 0.9922; %value from 30/05/2017 to 20/08/2018
else
    calibrationTool.tWindow = 0.9980; %value since 21/08/2018
end
