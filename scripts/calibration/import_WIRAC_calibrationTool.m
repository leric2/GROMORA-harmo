function calibrationTool = import_WIRAC_calibrationTool(calibrationTool)
%==========================================================================
% NAME      | import_GROMOS_calibrationTool.m
% TYPE      | Function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 12.2020
%           |
% ABSTRACT  | Complete the calibrationTool structure for GROMOS
%           |
% ARGUMENTS | INPUTS: 	1. calibrationTool: the default toolbox
%           |
%           | OUTPUTS: 	2. calibrationTool: the default toolbox completer for GROMOS
%           |
% COMMENTS  | External documentation for this structure is available on the
%           | git server of IAP
%           |
%==========================================================================

% Check that instrument name corresponds to this function
assert(strcmp(calibrationTool.instrumentName, 'WIRAC'),'Wrong instrument toolbox !')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.dataLocation='BERN';
calibrationTool.PI_NAME='Murk;Axel';
calibrationTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
calibrationTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
calibrationTool.PI_EMAIL='axel.murk@iap.unibe.ch';
calibrationTool.dataSource='MWR.O3_UBERN';

% Geolocation
calibrationTool.lon=7.44;
calibrationTool.lat=46.95;
calibrationTool.altitude=560;
calibrationTool.azimuthAngle=45;

calibrationTool.timeZone = 'Z';
calibrationTool.dateTime.TimeZone = calibrationTool.timeZone;

% Observation frequency
calibrationTool.observationFreq=1.4217504e11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrometer data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.numberOfSpectrometer=1;
calibrationTool.spectrometer='AC240';
calibrationTool.samplingRateFFTS=2000; % sampling rates in MHz 

calibrationTool.numberOfChannels=32768;

calibrationTool.IQProcessing = false;

% Local oscillators information
calibrationTool.fLO1=1.45875e11;
calibrationTool.fLO2=3.6e9;

% This one should correspond to the DC channel
calibrationTool.LOFreqTot=calibrationTool.fLO1-calibrationTool.fLO2-0.5e9;

%calibrationTool.DCChannel=16384;
calibrationTool.instrumentBandwidth=1e9;
%calibrationTool.LOFreq=1.45875e11;

% Known bad channels for the instrument
calibrationTool.badChannels=[16384 16385];

calibrationTool.numberOfAquisitionSpectraHot=60000;
calibrationTool.numberOfAquisitionSpectraAntenna=120000;
calibrationTool.numberOfAquisitionSpectraCold=120000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folder, Raw and log file data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path definition (for local computer only)
% calibrationTool.rawFileFolder=['/scratch/GROMOS_rawData/' calibrationTool.dateStr(1:4) '/' calibrationTool.dateStr(6:7) '/'];
% taken on the IAP lake, To Be mounted beforehand

% Valid properties for all instruments
calibrationTool.binaryDataExtension = '.bin';
calibrationTool.logFileDataExtension = '.txt';

calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';

calibrationTool.rawFileFolder=['/storage/lake/instrumentdata/wirac/' calibrationTool.dateStr(1:4) '/'];
%calibrationTool.extraFileFolder='/storage/atmosphere/instruments/gromos/level1/GROMORA/ExtraRawFiles/'; % no write permission on the IAP lake
calibrationTool.level1Folder='/home/esauvageat/Documents/WIRAC/Level1/';
%calibrationTool.level1Folder='/home/eric/Documents/PhD/GROSOM/Data/Level1/';

switch calibrationTool.dateStr 
    case '2022_02_28'
        calibrationTool.filename='wira-c_20220228_TC1';
    case '2022_03_01'
        calibrationTool.filename='wira-c_20220301';
    case '2022_03_02'
        calibrationTool.filename='wira-c_20220302_TC';
end
calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];

% Defining level1a filename to read (to be adapted for other users)
% calibrationTool.filenameLevel1a = [calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '_Tcold.nc'];

calibrationTool.checkLevel0=false;

% Log file
calibrationTool.delimiter_logfile = '\t';
calibrationTool.THotUnit='degreeC';

% Function for the harmonization of the log
calibrationTool.harmonize_log=@(calibrationtTool, log) harmonize_log_wirac(calibrationtTool, log);

calibrationTool.positionIndAsName = false;

calibrationTool.indiceCold=2;
calibrationTool.indiceAntenna=4;
calibrationTool.indiceHot=1;
calibrationTool.indiceTC =9;

calibrationTool.elevationAngleAntenna=140; %40
calibrationTool.elevationAngleCold=280; %-84
calibrationTool.elevationAngleHot=90;%

calibrationTool.elevationAngleTolerance=0;
calibrationTool.elevationAngleHotTol = 0;
calibrationTool.elevationAngleColdTol = 15;

calibrationTool.cycleDurationCold = 10;
calibrationTool.cycleDurationSky = 7;
calibrationTool.cycleDurationHot = 10;

calibrationTool.flipped_spectra=false;
calibrationTool.flipAroundChannel = 16384;
calibrationTool.flip_spectra=@(rawSpectra, calibrationTool) flip_spectra_gromos(rawSpectra, calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw File and Log
calibrationTool.numberOfTippingCurveExpected=48;
calibrationTool.toleranceTippingCurves=2;

calibrationTool.goodFlagLN2Above = 1;
calibrationTool.goodFlagLN2Below = 0;

% Considering the expected number of tipping curve:
calibrationTool.numberOfCyclesExpected=1500;
calibrationTool.toleranceNumberCycles=0.01*calibrationTool.numberOfCyclesExpected;
calibrationTool.tippingSize=27;

% Temperatures
calibrationTool.TNoiseCenterTh=650;
calibrationTool.TNoiseThresh=50;
calibrationTool.stdTNoiseThresh=8;

calibrationTool.THotTh=313.9;
calibrationTool.THotAbsThresh=2;

calibrationTool.hotTemperatureStdThreshold=0.05;

% Calibration
% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid
calibrationTool.stdAntAngleThresh = 0.5;
calibrationTool.adcOverloadThresh = 0;

calibrationTool.minNumberOfIndicePerCycle=12;
calibrationTool.threshNumRawSpectraHot=0.1*calibrationTool.numberOfChannels;
calibrationTool.threshNumRawSpectraCold=0.1*calibrationTool.numberOfChannels;
calibrationTool.threshNumRawSpectraAnt = 0.1*calibrationTool.numberOfChannels;

calibrationTool.maxProportionOfIndLN2LevelOutlier = 0.2;
calibrationTool.maxProportionOfIndLN2SensorOutlier = 0.2;

calibrationTool.frequencyBandAroundCenterTNoise = 200e6;

% Frequency lock flag
calibrationTool.maxProportionFreqLockError = 0.1;

% Max std dev of Gunn voltage
calibrationTool.maxStdV_Gun = 1e-1;

% Filters for flagging "bad channels"
calibrationTool.maxStdDevTbCal = 25; %TODO
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
calibrationTool.meteoFolder=['/storage/lake/instrumentdata/meteo/exwi/meteo/' calibrationTool.dateStr(1:4) '/'];
%calibrationTool.meteoFolder='/home/eric/Documents/PhD/GROSOM/Data/METEO_DATA/';

% Read meteo data
calibrationTool.read_meteo_data =@(calibrationTool) read_meteo_data_unibe(calibrationTool);

% Add meteo data to calibrated spectra
% calibrationTool.add_meteo_data = @(calibrationTool,correctedSpectra) add_meteo_data_unibe(calibrationTool,correctedSpectra);
calibrationTool.add_meteo_data = @(calibrationTool, meteoData, correctedSpectra) add_meteo_data_generic(calibrationTool, meteoData, correctedSpectra);

% Backup reading of MCH ground station data (ANETZ)
calibrationTool.meteoAnetzFolder = ['/storage/atmosphere/MeteoSchweiz/' calibrationTool.dateStr(1:4) '/'];
calibrationTool.anetzStnName = 'BER';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tipping curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO
calibrationTool.doTippingCurve = true;
%calibrationTool.TC.type = 'SkyLoads';
calibrationTool.TC.numberOfChannelsTropCorr = 500;
calibrationTool.TC.skipFraction = 0.05;
calibrationTool.TC.useWings = 'both';
calibrationTool.TC.deltaT = 10.4;
calibrationTool.TC.tauInitTC = 0.3;
calibrationTool.TC.maxIterTC = 500;
calibrationTool.TC.offsetTC = 5e-2;

% tipping curve
calibrationTool.run_tipping_curve = @(rawSpectra, log, calibrationTool) run_tipping_curve_generic(rawSpectra,log, calibrationTool);
calibrationTool.get_tipping_curve_data = @(rawSpectra, log, calibrationTool) get_tipping_curve_data_generic(rawSpectra,log, calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrections
calibrationTool.transmittanceWindow=0.99;

% Corrections
calibrationTool.troposphericCorrection.type = 'Ingold_v1_fit';
calibrationTool.troposphericCorrection.useWings = 'both';
calibrationTool.troposphericCorrection.numberOfChannelsTropCorr = 500;
calibrationTool.troposphericCorrection.skipFraction = 0.05;
calibrationTool.troposphericCorrection.deltaT = 10.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level0 -> Level1a functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the functions that will be used for processing this
% calibration

% Reading routine to use for the raw data
calibrationTool.read_level0=@(calibrationTool, readRawFile) read_level0_wirac(calibrationTool, readRawFile);

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
calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,N) plot_spectra_generic_nodrift(calibrationTool,drift,meteoData, calibratedSpectra,N);

% Function for quality check of the calibrated spectra
calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_wirac(log,calibrationTool,calibratedSpectra);

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
% binary name extension

end
