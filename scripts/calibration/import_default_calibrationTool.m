function calibrationTool = import_default_calibrationTool(instrumentName,dateStr)
%==========================================================================
% NAME          | import_default_calibrationTool(instrumentName)
% TYPE          | Function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: - instrumentName: name of the instrument as a
%               | string (ex: 'GROMOS')
%               |
%               | OUTPUTS: - calibrationTool: the default toolbox for
%               | launching a retrieval for this instrument.
%               | 
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
calibrationTool = struct();

% dateStr
calibrationTool.dateStr=dateStr;

calibrationTool.Year = str2double(dateStr(1:4));
calibrationTool.Month = str2double(dateStr(6:7));
calibrationTool.Day = str2double(dateStr(9:10));

calibrationTool.dateTime = datetime(dateStr,'InputFormat','yyyy_MM_dd');

% Name of the instrument
calibrationTool.instrumentName=instrumentName;
calibrationTool.dateStr=dateStr;

% Valid properties for all instruments
calibrationTool.binaryDataExtension = '.bin';
calibrationTool.logFileDataExtension = '.txt';

calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';

calibrationTool.positionIndAsName = false;

% Physical constant
calibrationTool.lightSpeed=299792458; % [m/s] 
calibrationTool.h=6.62606876e-34; % [J/s]
calibrationTool.kb=1.38065e-23;    % [J/K]

calibrationTool.zeroDegInKelvin = 273.15;

calibrationTool.backgroundMWTb = 2.7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing default parameters for each instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROMOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch instrumentName
    case 'GROMOS'
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
        
        % Observation frequency
        calibrationTool.observationFreq=1.4217504e11;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        calibrationTool.numberOfSpectrometer=1;
        calibrationTool.spectrometer='AC240';
        calibrationTool.samplingRateFFTS=2000;
       
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
        % calibrationTool.rawFileFolder=['/scratch/GROMOS_rawData/' dateStr(1:4) '/' dateStr(6:7) '/'];
        % taken on the IAP lake, To Be mounted beforehand
        calibrationTool.rawFileFolder=['/mnt/instrumentdata/gromos/FFTS/' dateStr(1:4) '/'];
        %calibrationTool.rawFileFolder=['/home/eric/Documents/PhD/GROSOM/rawData/'];
        calibrationTool.extraFileFolder='/scratch/GROSOM/ExtraRawFiles/'; % no write permission on the IAP lake
        calibrationTool.level1Folder='/scratch/GROSOM/Level1/GROMOS/';
        %calibrationTool.level1Folder='/home/eric/Documents/PhD/GROSOM/Level1/';
        
        calibrationTool.filename=[calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
        calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];
        
        calibrationTool.checkLevel0=true;
        
        % Log file
        calibrationTool.delimiter_logfile = '\t';
        calibrationTool.THotUnit='degreeC';
        
        % Function for the harmonization of the log
        calibrationTool.harmonize_log=@(calibrationtTool, log) harmonize_log_gromos(calibrationtTool, log);  
        
        calibrationTool.indiceCold=0;
        calibrationTool.indiceAntenna=1;
        calibrationTool.indiceHot=2;
        
        calibrationTool.elevationAngleAntenna=40;
        calibrationTool.elevationAngleCold=-84;
        calibrationTool.elevationAngleHot=160;
        
        calibrationTool.elevationAngleTolerance=2;
        calibrationTool.elevationAngleHotTol = 0;
        calibrationTool.elevationAngleColdTol = 0;
        
        calibrationTool.cycleDurationCold = 10;
        calibrationTool.cycleDurationSky = 7;
        calibrationTool.cycleDurationHot = 10;
        
        calibrationTool.flipped_spectra=true;
        calibrationTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);   
        
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
        calibrationTool.TSysCenterTh=2750;
        calibrationTool.TSysThresh=300;
        calibrationTool.stdTSysThresh=8;
        
        calibrationTool.THotTh=313.9;
        calibrationTool.THotAbsThresh=2;
        
        calibrationTool.hotTemperatureStdThreshold=0.05;
        
        % Calibration
        % minimum number of indices (h-a-c) we want in a calibration cycle for it
        % to be valid
        calibrationTool.stdAntAngleThresh = 0.5;
        
        calibrationTool.minNumberOfIndicePerCycle=12;
        calibrationTool.threshNumRawSpectraHot=0.1*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraCold=0.1*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraAnt = 0.1*calibrationTool.numberOfChannels;
        
        calibrationTool.maxProportionOfIndLN2LevelOutlier = 0.2;
        calibrationTool.maxProportionOfIndLN2SensorOutlier = 0.2;
        
        calibrationTool.frequencyBandAroundCenterTSys = 200e6;
        
        % Filters for flagging "bad channels"
        calibrationTool.maxStdDevTb = nan; %TODO
        
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
        calibrationTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        %calibrationTool.meteoFolder='/home/eric/Documents/PhD/GROSOM/METEO_DATA/';
        
        % Read meteo data
        calibrationTool.read_meteo_data =@(calibrationTool) read_meteo_data_unibe(calibrationTool);
        
        % Add meteo data to calibrated spectra
        % calibrationTool.add_meteo_data = @(calibrationTool,correctedSpectra) add_meteo_data_unibe(calibrationTool,correctedSpectra);
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
        calibrationTool.tWindow=0.99;
        
        % Corrections
        calibrationTool.troposphericCorrection.type = 'Ingold_v1';
        calibrationTool.troposphericCorrection.useWings = 'both';
        calibrationTool.troposphericCorrection.numberOfChannelsTropCorr = 50;
        calibrationTool.troposphericCorrection.skipFraction = 0.05;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Level0 -> Level1a functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Selecting the functions that will be used for processing this
        % calibration
        
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
        calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N) plot_spectra_generic(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N);
        
        % Function for quality check of the calibrated spectra
        calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_generic(log,calibrationTool,calibratedSpectra);
        
        % Function saving the calibrated spectra into netCDF file
        calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(calibrationTool,log,calibratedSpectra,warningLevel0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Level1a -> Level1b functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Function reading the daily calibrated spectra from netCDF file
        calibrationTool.read_level1a = @(calibrationTool) read_level1a_daily(calibrationTool);
        
        % Check of the channels quality on the calibrated spectra:
        calibrationTool.checking_channel_quality= @(calibratedSpectra,calibrationTool,filterN) checking_channel_quality_gromos(calibratedSpectra,calibrationTool,filterN);
        
        % Integration of level1a data
        calibrationTool.integrate_calibrated_spectra= @(calibrationTool,calibratedSpectra) integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra);
        
        % Function for plotting the integrated spectra (when hourly)
        calibrationTool.plot_integrated_spectra = @(calibrationTool,rawSpectra,lowerLim,upperLim) plot_integrated_spectra_generic(calibrationTool,rawSpectra,lowerLim,upperLim);
        
        calibrationTool.tropospheric_correction = @(integration,calibrationTool,TtropCorr) tropospheric_correction_generic(integration,calibrationTool,TtropCorr);
        
        % Window correction for the calibrated spectra
        calibrationTool.window_correction= @(calibrationTool,level1b) window_correction_generic(calibrationTool,level1b);
        
        % Check of the integrated spectra
        calibrationTool.check_integrated = @(calibrationTool,integratedSpectra) check_integrated_generic(calibrationTool,integratedSpectra);
        
        % Function saving the calibrated spectra into netCDF file
        calibrationTool.save_level1b=@(calibrationTool,level1b) save_level1b_daily(calibrationTool,level1b);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOMORA    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SOMORA'
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
        calibrationTool.rawFileFolder=['/scratch/SOMORA_rawData/2019/' dateStr(6:7) '/'];
        
        %calibrationTool.rawFileFolder=['/home/eric/Documents/PhD/GROSOM/rawData/'];
        %calibrationTool.level1Folder='/home/esauvageat/Documents/GROSOM/Analysis/Level1/SOMORA/';
        calibrationTool.extraFileFolder='/scratch/GROSOM/ExtraRawFiles/'; % no write permission on the IAP lake
        %calibrationTool.level1Folder='/home/eric/Documents/PhD/GROSOM/Level1/';
        calibrationTool.level1Folder='/scratch/GROSOM/Level1/SOMORA/';
        calibrationTool.filename=[calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
        calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];
        
        % Log
        calibrationTool.delimiter_logfile = '\t';
        calibrationTool.harmonize_log=@(calibrationTool, log) harmonize_log_somora(calibrationTool, log);
        calibrationTool.THotUnit='K';
        
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
        calibrationTool.TSysCenterTh=2750;
        calibrationTool.TSysThresh=300;
        calibrationTool.stdTSysThresh=8;
        
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
        
        calibrationTool.frequencyBandAroundCenterTSys = 200e6;
        
        calibrationTool.maxProportionOfIndLN2LevelOutlier = 0.3;
        calibrationTool.maxProportionOfIndLN2SensorOutlier = 0.3;
        
        % Filters for flagging "bad channels"
        calibrationTool.maxStdDevTb = 10;
        
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
        calibrationTool.meteoFolder=['/scratch/GROSOM/MeteoDataSOMORA/METEO_DATA_' dateStr(1:4) '/'];
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
        calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N) plot_spectra_generic(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N);
        
        % Function for quality check of the calibrated spectra
        calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_generic(log,calibrationTool,calibratedSpectra);
        
        % Function saving the calibrated spectra into netCDF file
        calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_level1a_daily(calibrationTool,log,calibratedSpectra,warningLevel0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Level1a -> Level1b functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Function reading the daily calibrated spectra from netCDF file
        calibrationTool.read_level1a = @(calibrationTool) read_level1a_daily(calibrationTool);
        
        % Check of the channels quality on the calibrated spectra:
        calibrationTool.checking_channel_quality= @(calibratedSpectra,calibrationTool,filterN) checking_channel_quality_gromos(calibratedSpectra,calibrationTool,filterN);
        
        % Integration of level1a data
        calibrationTool.integrate_calibrated_spectra= @(calibrationTool,calibratedSpectra) integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra);
        
        % Function for plotting the integrated spectra (when hourly)
        calibrationTool.plot_integrated_spectra = @(calibrationTool,rawSpectra,lowerLim,upperLim) plot_integrated_spectra_generic(calibrationTool,rawSpectra,lowerLim,upperLim);
        
        calibrationTool.tropospheric_correction = @(integration,calibrationTool,TtropCorr) tropospheric_correction_generic(integration,calibrationTool,TtropCorr);
        
        % Window correction for the calibrated spectra
        calibrationTool.window_correction= @(calibrationTool,level1b) window_correction_generic(calibrationTool,level1b);
        
        % Check of the integrated spectra
        calibrationTool.check_integrated = @(calibrationTool,integratedSpectra) check_integrated_generic(calibrationTool,integratedSpectra);
        
        % Function saving the calibrated spectra into netCDF file
        calibrationTool.save_level1b=@(calibrationTool,level1b) save_level1b_daily(calibrationTool,level1b);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOPI       
    case 'mopi5'
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
        calibrationTool.rawFileFolder=['/mnt/instrumentdata/mopi5/' dateStr(1:4) '/'];
        calibrationTool.level1Folder='/scratch/MOPI5/Level1/';
        calibrationTool.extraFileFolder='/scratch/MOPI5/ExtraRawFiles/'; % no write permission on the IAP lake
        
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
        calibrationTool.stdAntAngleThresh = 0.5;
        
        calibrationTool.minNumberOfIndicePerCycle=15;

        calibrationTool.threshNumRawSpectraHot=0.1*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraCold=0.1*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraAnt = 0.1*calibrationTool.numberOfChannels;
        
        calibrationTool.maxProportionOfIndLN2LevelOutlier = 0.3;
        calibrationTool.maxProportionOfIndLN2SensorOutlier = 0.3;
        
        % Filters for flagging "bad channels"
        % On 10 minutes spectra
        calibrationTool.maxStdDevTb = nan; %TODO
        calibrationTool.maxStdDevTbInt = 4;
        
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
        calibrationTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        %calibrationTool.meteoFolder='/home/eric/Documents/PhD/METEO_DATA/';
    
        % Function specific to this instrument
        % meteo Data
        %calibrationTool.get_meteo_data = @(calibrationTool,correctedSpectra) get_meteo_data_payerne(calibrationTool,correctedSpectra);
        calibrationTool.read_meteo_data =@(calibrationTool) read_meteo_data_unibe(calibrationTool);
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
        calibrationTool.tWindow=0.99;
        
        calibrationTool.troposphericCorrection.type = 'Ingold_v1_fit';
        calibrationTool.troposphericCorrection.useWings = 'both';
        calibrationTool.troposphericCorrection.numberOfChannelsTropCorr = 100;
        calibrationTool.troposphericCorrection.skipFraction = 0.05;

        
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
        calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N) plot_spectra_generic(calibrationTool,drift,meteoData, calibratedSpectra,lowerLim,upperLim,N);

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
        calibrationTool.checking_channel_quality= @(calibratedSpectra,calibrationTool,filterN) checking_channel_quality_gromos(calibratedSpectra,calibrationTool,filterN);
        
        % Integration of level1a data
        calibrationTool.integrate_calibrated_spectra= @(calibrationTool,calibratedSpectra) integrate_calibrated_spectra_generic(calibrationTool,calibratedSpectra);
        
        % Function for plotting the integrated spectra (when hourly)
        calibrationTool.plot_integrated_spectra = @(calibrationTool,rawSpectra,lowerLim,upperLim) plot_integrated_spectra_generic(calibrationTool,rawSpectra,lowerLim,upperLim);
        
        calibrationTool.tropospheric_correction = @(integration,calibrationTool,TtropCorr) tropospheric_correction_generic(integration,calibrationTool,TtropCorr);
        
        % Window correction for the calibrated spectra
        calibrationTool.window_correction= @(calibrationTool,level1b) window_correction_generic(calibrationTool,level1b);
        
        % Check of the integrated spectra
        calibrationTool.check_integrated = @(calibrationTool,integratedSpectra) check_integrated_generic(calibrationTool,integratedSpectra);
        
        % Function saving the calibrated spectra into netCDF file
        calibrationTool.save_level1b=@(calibrationTool,level1b) save_level1b_daily(calibrationTool,level1b);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIAWARA-C
    case 'MIAWARA-C'
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
    
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter varying with time for the instruments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conversion of datestr to datenum:
timeNumber=datenum(str2num(dateStr(1:4)),str2num(dateStr(6:7)),str2num(dateStr(9:10)));

switch instrumentName
    case 'GROMOS'
        % binary name extension 
        if timeNumber < datenum(2010,05,06)
            calibrationTool.binaryDataExtension = '.dat';
        end
        if timeNumber < datenum(2010,03,10)
            calibrationTool.positionIndAsName = true;
            calibrationTool.nameColdIndice = 'Cold';
            calibrationTool.nameHotIndice = 'Hot';
            calibrationTool.nameAntennaIndice = 'Antenna';
            calibrationTool.otherName = 'Reference';
        end
        
        if timeNumber < datenum(2012,01,01) %TOCHECK
            calibrationTool.goodFlagLN2Above = 1;
            calibrationTool.goodFlagLN2Below = 1;
            
            calibrationTool.stdTSysThresh = 15;
        end
        % window transmission 
        if timeNumber < datenum(2018,11,12)
            calibrationTool.tWindow = 0.99; % has been changed at that time but no idea of the values ??????
        else
            calibrationTool.tWindow = 0.99; 
        end
        
        % Elevation angle of the cold load
        if  timeNumber < datenum(2012,04,26) 
            % Before April 2012, same angle measured for the 3 position
            calibrationTool.elevationAngleAntenna=39.85;
            calibrationTool.elevationAngleCold=39.85;
            calibrationTool.elevationAngleHot=39.85;
            calibrationTool.elevationAngleTolerance= 2;
            calibrationTool.elevationAngleHotTol = 1;
            calibrationTool.elevationAngleColdTol = 1;
        elseif (timeNumber>= datenum(2012,04,26) && timeNumber< datenum(2013,01,29))
            calibrationTool.elevationAngleCold=-85;
        elseif (timeNumber>= datenum(2013,01,29) && timeNumber<datenum(2014,09,19)) 
             calibrationTool.elevationAngleCold=-84;
        elseif (timeNumber>= datenum(2014,09,19) && timeNumber<datenum(2019,2,11)) 
                calibrationTool.elevationAngleCold=-85;
        elseif (timeNumber>= datenum(2019,02,12) && timeNumber<datenum(2019,3,12)) 
                calibrationTool.elevationAngleCold=-89;
        elseif timeNumber >= datenum(2019,3,12)
            calibrationTool.elevationAngleCold=-84;
        end
        
         if (timeNumber>= datenum(2019,01,14) && timeNumber<datenum(2019,03,12)) 
             calibrationTool.elevationAngleBias = 5;
         end
        
    case 'SOMORA'
         if timeNumber>datenum(2016,01,01) && timeNumber<datenum(2017,04,01)
             calibrationTool.THotTh = 297.4;
         end
        
        if timeNumber<datenum(2019,01,01)
            calibrationTool.rawFileFolder=['/media/esauvageat/INTENSO/RAW_DATA/' dateStr(1:4) '/' dateStr(6:7) '/'];
            calibrationTool.filename=[calibrationTool.instrumentName,'09_', calibrationTool.dateStr];
            calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];
        end
        % window transmission 
        if timeNumber<datenum(2007,07,03) %attn anciennement  733957 cad 3 juillet 2009 mais FAUX
            calibrationTool.tWindow = 0.987; %value from 2002 to 3/7/2007 (att. window has been changed 8/2/2006 but new t not measured) nb: ema 9/11/2007
        elseif (timeNumber>= datenum(2007,07,03) && timeNumber<datenum(2010,6,26)) %ie 734315
            calibrationTool.tWindow = 0.979; %value from 3/7/2007 to 25/6/2010   nb: ema 9/11/2007; calc sur la base de R-J
        elseif (timeNumber>= datenum(2010,06,25) && timeNumber<datenum(2012,03,02)) %ie 734931 2/3/2012
            calibrationTool.tWindow = 0.985; %value from 26/6/2010 to 2/3/2012 %attn anciennement 0.976 ???
        elseif (timeNumber>= datenum(2012,03,02) && timeNumber<datenum(2014,06,16)) %ie 16/06/2014
            calibrationTool.tWindow = 0.999; %value from 3/3/2012 to 15/06/2014
        elseif (timeNumber>= datenum(2014,06,16) && timeNumber<datenum(2015,09,08))
            calibrationTool.tWindow = 1.0105; %value from 16/6/2014 to 07/09/2015
        elseif (timeNumber>= datenum(2015,09,09) && timeNumber<datenum(2017,05,30))
            calibrationTool.tWindow = 1.009; %value from 8/09/2015 to 29/05/2017
        elseif (timeNumber>= datenum(2017,05,31) && timeNumber<datenum(2018,08,21))
            calibrationTool.tWindow = 0.9922; %value from 30/05/2017 to 20/08/2018
        else
            calibrationTool.tWindow = 0.9980; %value since 21/08/2018
        end
        
    case 'mopi5'
        if timeNumber<datenum(2019,02,12)
            calibrationTool.elevationAngleHot=90;
        end
    case 'MIAWARA-C'
        % Meteo
        YYYY = calibrationTool.dateStr(1:4);
        MM   = calibrationTool.dateStr(6:7);
        DD   = calibrationTool.dateStr(9:10);
        
        dat  = datenum(dateStr);
        if dat > datenum('2015-09-01') && dat < datenum(2016,11,25)
            calibrationTool.get_meteo_data = @(calibrationTool) get_meteo_data_grc(calibrationTool);
            %calibrationTool.meteoFile = ['/data/miradara3/instrumentdata/gromosc/' YYYY '/GROMOS-C*' YYYY '_' MM '_' DD];
            calibrationTool.meteoFile = [calibrationTool.rawFileFolder 'GROMOS-C*' YYYY '_' MM '_' DD];
        elseif dat > datenum('2017-04-06 23:00')
            calibrationTool.get_meteo_data = @(calibrationTool) get_meteo_data_miac_ownfile(calibrationTool);
            calibrationTool.meteoFile = [ calibrationTool.rawFileFolder 'MIAWARA-C_meteo_' YYYY MM DD];

        else
            calibrationTool.get_meteo_data = @(log) get_meteo_data_miac(log);
        end
end
