function calibrationTool = import_GROMOS_FB_calibrationTool(calibrationTool)
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
assert(strcmp(calibrationTool.instrumentName, 'GROMOS'),'Wrong instrument toolbox !')

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
calibrationTool.spectrometer='FB';
calibrationTool.samplingRateFFTS=2000; % sampling rates in MHz 

calibrationTool.numberOfChannels=48;

calibrationTool.IQProcessing = false;

% Frequencies theory and measurement
calibrationTool.channel_width_theory = 1e6*[100 100 100 100 100 100 30 20 20 10 5 5 2 2 2 1 1 0.5 0.5 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.5 0.5 1 1 2 2 2 5 5 10 20 20 30 100 100 100 100 100 100 100];
calibrationTool.channel_width = 1e6*[114.667784	118.600345	109.055542	146.597255	164.5163	18.906618	18.189048	21.780333	10.595859	5.630094	6.079148	2.830508	2.815202	3.103099	1.794136	1.39955	0.826034	0.770904	0.242118	0.230111	0.225406	0.231984	0.243707	0.23766	0.231579	0.235805	0.233243	0.646687	0.642863	1.560209	1.261859	3.162788	2.869584	3.370188	6.127845	7.779886	11.394174	22.580517	23.252718	31.629541	178.777324	170.339139	180.648311	139.313079	159.515409	180.338123	151.604634];

calibrationTool.peak_voltage = [4.473	4.056	3.961	2.987	2.454	0.136	0.21	0.6	0.821	0.746	0.654	1.044	1.05	0.96	1.846	2.47	3.723	4.587	2.776	2.94	2.991	2.89	2.687	2.799	2.804	2.788	2.956	5.532	4.748	2.248	2.81	1.259	1.161	1	0.58	0.527	1.15	0.58	0.607	0.501	1.958	1.952	2.161	2.699	2.239	2.088	2.719];

calibrationTool.intermediate_freq_theory = [3700000000.0 4265960000.0 4169960000.0 4067960000.0 3978960000.0 3851960000.0 3779960000.0 3761560000.0 3739960000.0 3726560000.0 3716960000.0 3713060000.0 3708960000.0 3706860000.0 3704760000.0 3703460000.0 3702560000.0 3701360000.0 3700960000.0 3700760000.0 3700560000.0 3700360000.0 3700160000.0 3700000000.0 3699780000.0 3699560000.0 3699360000.0 3699160000.0 3698760000.0 3698260000.0 3697660000.0 3696660000.0 3695160000.0 3692760000.0 3691560000.0 3687260000.0 3684760000.0 3674860000.0 3660960000.0 3640360000.0 3616460000.0 3559960000.0 3433960000.0 3329960000.0 3230960000.0 3144960000.0 3650000000.0 3750000000.0]-3700e6;
calibrationTool.intermediate_freq = [-9999, 569.699054	477.434267	372.997632	280.70649	180.975139	91.163342	57.83471	38.991837	26.74233	16.691558	12.855815	8.802827	6.662681	4.756655	3.406094	2.607902	1.305915	1.055343	0.729692	0.541362	0.350396	0.150986	-0.057858	-0.25199	-0.443818	-0.653617	-0.840308	-1.301629	-1.687485	-2.346339	-3.275542	-4.874344	-7.358491	-8.439596	-13.013544	-15.636809	-25.483625	-39.911915	-59.739096	-85.378971	-164.135651	-270.444599	-367.113031	-471.926288	-550.152046	-40.215867	74.744203];

% Local oscillators information
calibrationTool.fLO1=145.87504e9; %1.45875e11;
calibrationTool.fLO2=3.7e9;

% This one should correspond to the DC channel
%calibrationTool.LOFreqTot=calibrationTool.fLO1-calibrationTool.fLO2-0.5e9;
calibrationTool.LOFreqTot=calibrationTool.fLO1-calibrationTool.fLO2;

%calibrationTool.DCChannel=16384;
calibrationTool.instrumentBandwidth=1e9;
%calibrationTool.LOFreq=1.45875e11;

% Known bad channels for the instrument
calibrationTool.badChannels=[1 41 47 48];

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

calibrationTool.rawFileFolder=['/storage/lake/instrumentdata/gromos/FFTS/' calibrationTool.dateStr(1:4) '/'];
%calibrationTool.rawFileFolder=['/home/eric/Documents/PhD/GROSOM/Data/rawData/'];
calibrationTool.extraFileFolder='/storage/tub/instruments/gromos/level1/GROMORA/ExtraRawFiles/'; % no write permission on the IAP lake
calibrationTool.level1Folder=['/storage/tub/instruments/gromos/level1/GROMORA/v2/' calibrationTool.dateStr(1:4) '/'];
%calibrationTool.level1Folder='/home/eric/Documents/PhD/GROSOM/Data/Level1/';

calibrationTool.filename=[calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];

% Defining level1a filename to read (to be adapted for other users)
%calibrationTool.filenameLevel1a = [calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '_Tcold.nc'];

calibrationTool.checkLevel0=true;

% Log file
calibrationTool.delimiter_logfile = '\t';
calibrationTool.THotUnit='degreeC';

% Function for the harmonization of the log
calibrationTool.harmonize_log=@(calibrationtTool, log) harmonize_log_gromos(calibrationtTool, log);

calibrationTool.positionIndAsName = false;

calibrationTool.indiceCold=0;
calibrationTool.indiceAntenna=1;
calibrationTool.indiceHot=2;
calibrationTool.indiceTC = 5;

calibrationTool.elevationAngleAntenna=40;
calibrationTool.elevationAngleCold=-84;
calibrationTool.elevationAngleHot=160;

calibrationTool.elevationAngleTolerance=2;
calibrationTool.elevationAngleHotTol = 0;
calibrationTool.elevationAngleColdTol = 0;

calibrationTool.cycleDurationCold = 1;
calibrationTool.cycleDurationSky = 1;
calibrationTool.cycleDurationHot = 1;

calibrationTool.flipped_spectra=true;
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
calibrationTool.TNoiseCenterTh=2750;
calibrationTool.TNoiseThresh=300;
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

calibrationTool.filterTypeChannelQualityCal = 5;
calibrationTool.filterTypeChannelQualityInt = 5;

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
calibrationTool.add_meteo_data = @(calibrationTool, meteoData, correctedSpectra) add_meteo_data_GROMOS_FB(calibrationTool, meteoData, correctedSpectra);
calibrationTool.meteoTimeExtension = 10;

% Backup reading of MCH ground station data (ANETZ)
calibrationTool.meteoAnetzFolder = ['/storage/tub/MeteoSchweiz/' calibrationTool.dateStr(1:4) '/'];
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
calibrationTool.troposphericCorrection.numberOfChannelsTropCorr = 0;
calibrationTool.troposphericCorrection.skipFraction = 0.05;
calibrationTool.troposphericCorrection.deltaT = 10.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level0 -> Level1a functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the functions that will be used for processing this
% calibration

% Reading routine to use for the raw data
calibrationTool.read_level0=@(calibrationTool, readRawFile) read_level0_old_FFTS_GROMOS(calibrationTool, readRawFile);

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
calibrationTool.read_level1a = @(calibrationTool) read_level1_FB_GROMORA(calibrationTool);

% Check of the channels quality on the calibrated spectra:
calibrationTool.check_channel_quality= @(calibratedSpectra,calibrationTool,filterN) check_channel_quality_generic(calibratedSpectra,calibrationTool,filterN);

% Integration of level1a data
calibrationTool.integrate_calibrated_spectra= @(calibrationTool,calibratedSpectra) integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra);

% Function for plotting the integrated spectra (when hourly)
calibrationTool.plot_integrated_spectra = @(calibrationTool,rawSpectra) plot_integrated_spectra_FB(calibrationTool,rawSpectra);

calibrationTool.tropospheric_correction = @(integration,calibrationTool) tropospheric_correction_FB(integration,calibrationTool);

% Window correction for the calibrated spectra
calibrationTool.window_correction= @(calibrationTool,level1b) window_correction_generic(calibrationTool,level1b);

% Check of the integrated spectra
calibrationTool.check_integrated = @(calibrationTool,integratedSpectra) check_integrated_FB(calibrationTool,integratedSpectra);

% Function saving the calibrated spectra into netCDF file
calibrationTool.save_level1b=@(calibrationTool,level1b) save_level1b_FB(calibrationTool,level1b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter varying with time for the instruments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% binary name extension
if calibrationTool.timeNumber < datenum(2010,05,06)
    calibrationTool.binaryDataExtension = '.dat';
end
if calibrationTool.timeNumber < datenum(2010,03,10)
    calibrationTool.positionIndAsName = true;
    calibrationTool.nameColdIndice = 'Cold';
    calibrationTool.nameHotIndice = 'Hot';
    calibrationTool.nameAntennaIndice = 'Antenna';
    calibrationTool.otherName = 'Reference';

    % Not working for now before 2010-03-10
    %calibrationTool.doTippingCurve = false;
end

% Noise Temperature
if calibrationTool.timeNumber <= datenum(2010,03,16)
     calibrationTool.TNoiseCenterTh = 2800;
    calibrationTool.stdTNoiseThresh = 15;
elseif (calibrationTool.timeNumber > datenum(2010,03,16)) && (calibrationTool.timeNumber <= datenum(2013,11,07))
    % Change in TNoise at that date
    calibrationTool.TNoiseCenterTh = 2550;
    calibrationTool.stdTNoiseThresh = 15;
elseif calibrationTool.timeNumber > datenum(2021,12,16)
    calibrationTool.TNoiseCenterTh = 3600;
    calibrationTool.stdTNoiseThresh = 15;
end

if calibrationTool.timeNumber <= datenum(2010,08,18)
    calibrationTool.goodFlagLN2Above = 1;
    calibrationTool.goodFlagLN2Below = 1;
elseif (calibrationTool.timeNumber > datenum(2010,08,18) && calibrationTool.timeNumber < datenum(2012,07,23)) %TOCHECK
    calibrationTool.goodFlagLN2Above = 1;
    calibrationTool.goodFlagLN2Below = 0;
    
    
elseif calibrationTool.timeNumber > datenum(2020,06,19)
    calibrationTool.goodFlagLN2Above = 0;
    calibrationTool.goodFlagLN2Below = 0;
end
% window transmission
if calibrationTool.timeNumber < datenum(2018,11,12)
    calibrationTool.transmittanceWindow = 0.99; % has been changed at that time but no idea of the values ??????
else
    calibrationTool.transmittanceWindow = 0.99;
end

% Elevation angle of the cold load
if  calibrationTool.timeNumber < datenum(2012,04,26)
    % Before 26.04.12, no drift quantities because the order of the
    % observation was not the same !
    calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,N) plot_spectra_generic_nodrift(calibrationTool,drift,meteoData, calibratedSpectra,N);
    calibrationTool.adcOverloadThresh = 10;

    % Before April 2012, same angle measured for the 3 position
    calibrationTool.elevationAngleAntenna=39.85;
    calibrationTool.elevationAngleCold=39.85;
    calibrationTool.elevationAngleHot=39.85;
    calibrationTool.elevationAngleTolerance= 2;
    calibrationTool.elevationAngleHotTol = 1;
    calibrationTool.elevationAngleColdTol = 1;
elseif (calibrationTool.timeNumber>= datenum(2012,04,26) && calibrationTool.timeNumber< datenum(2013,01,29))
    calibrationTool.elevationAngleCold=-85;
elseif (calibrationTool.timeNumber>= datenum(2013,01,29) && calibrationTool.timeNumber<datenum(2014,09,19))
    calibrationTool.elevationAngleCold=-84;
elseif (calibrationTool.timeNumber>= datenum(2014,09,19) && calibrationTool.timeNumber<datenum(2019,2,11))
    calibrationTool.elevationAngleCold=-85;
elseif (calibrationTool.timeNumber>= datenum(2019,02,12) && calibrationTool.timeNumber<datenum(2019,3,12))
    calibrationTool.elevationAngleCold=-89;
elseif calibrationTool.timeNumber >= datenum(2019,3,12)
    calibrationTool.elevationAngleCold=-84;
end

if (calibrationTool.timeNumber>= datenum(2019,01,14) && calibrationTool.timeNumber<datenum(2019,03,12))
    calibrationTool.elevationAngleBias = 5;
end

% TC

if (calibrationTool.timeNumber>= datenum(2010,01,01) && calibrationTool.timeNumber<datenum(2012,04,26))
    calibrationTool.indiceTC = 1;
end
if calibrationTool.timeNumber < datenum(2009,07,01)
    calibrationTool.read_level1a = @(calibrationTool) read_level1_FB_GROMORA(calibrationTool);
end
end


