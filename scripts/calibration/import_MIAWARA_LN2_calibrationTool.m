function calibrationTool = import_MIAWARA_LN2_calibrationTool(calibrationTool)
%==========================================================================
% NAME      | import_MIAWARAC_calibrationTool.m
% TYPE      | Function
% AUTHOR(S) | Eric Sauvageat; modified by A. Bell
% CREATION  | 01.2023
%           |
% ABSTRACT  | Complete the calibrationTool structure for MIAWARAC
%           |
% ARGUMENTS | INPUTS: 	1. calibrationTool: the default toolbox
%           |
%           | OUTPUTS: 	2. calibrationTool: the default toolbox completer for MIAWARA
%           |
% COMMENTS  | External documentation for this structure is available on the
%           | git server of IAP
%           |
%==========================================================================

%%%%%%%%%%%%%%%%%%Regularly changed fields%%%%%%%%%%%%%%%%%%%%%
suffix = '_LN2';
%suffix = '_manualmode';

calibrationTool.filename=[calibrationTool.instrumentName,'_', calibrationTool.dateStr(1:4) '_' calibrationTool.dateStr(6:7) '_' calibrationTool.dateStr(9:10) suffix];
calibrationTool.custom_mode = false;

% Check that instrument name corresponds to this function
assert(strcmp(calibrationTool.instrumentName, 'MIAWARA'),'Wrong instrument toolbox !')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta data
calibrationTool.dataLocation='Zimmerwald, Bern, CH';
calibrationTool.PI_NAME='Murk;Axel';
calibrationTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
calibrationTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
calibrationTool.PI_EMAIL='axel.murk@iap.unibe.ch';
calibrationTool.dataSource='MWR.H2O_UBERN002';

% Geolocation
calibrationTool.lon=7.465;
calibrationTool.lat= 46.877 ;
calibrationTool.altitude=907;
calibrationTool.azimuthAngle=0;
calibrationTool.flipAroundChannel = 0;

calibrationTool.timeZone = 'Z';
calibrationTool.dateTime.TimeZone = calibrationTool.timeZone;

% Observation frequency
calibrationTool.observationFreq=22.235;
calibrationTool.calibrationTime=1*60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrometer data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO - Add second spectrometer in proper way for l1a output

%calibrationTool.samplingRateFFTS=2000; % sampling rates in MHz 
if datenum(calibrationTool.dateStr) < datenum(2022,10,17)
    disp('date is before USRP installed')
    calibrationTool.numberOfSpectrometer=1;
    calibrationTool.spectrometer='AC240';
    calibrationTool.numberOfChannels = 16384;
    calibrationTool.IQProcessing = false;
    calibrationTool.badChannels = horzcat([1:64],[calibrationTool.numberOfChannels-64:calibrationTool.numberOfChannels]);
    calibrationTool.IQProcessing = false;
    cal_upgradeDONE = false;

else
    calibrationTool.numberOfSpectrometer=1;
    calibrationTool.numberOfChannels = 49152;
    calibrationTool.IQProcessing = false;
    cal_upgradeDONE = true;
    calibrationTool.spectrometer=['AC240' 'USRP1' 'USRP2'];
end

calibrationTool.instrumentBandwidth=1e9;

% Known bad channels for the instrument
calibrationTool.badChannels=[8192 8193];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folder, Raw and log file data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path definition (for local computer only)
% calibrationTool.rawFileFolder=['/scratch/GROMOS_rawData/' calibrationTool.dateStr(1:4) '/' calibrationTool.dateStr(6:7) '/'];
% taken on the IAP lake, To Be mounted beforehand

%valid for all instruments
calibrationTool.binaryDataExtension = '.bin';
calibrationTool.logFileDataExtension = '.txt';

calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';
calibrationTool.positionIndAsName = false;

%Paths
calibrationTool.extraFileFolder='/export/data/miawara/ExtraRawFiles/'; % no write permission on the IAP lake
calibrationTool.rawFileFolder=['/mnt/lake/instrumentdata/miawara/' calibrationTool.dateStr(1:4) '/'];%['/storage/lake/instrumentdata/miawarac/' dateStr(1:4) '/'];%['/mnt/instrumentdata/miawarac/' dateStr(1:4) '/'];%AB
calibrationTool.level1Folder = '/home/alistair/export/data/miawara/L2_nc_Retrievals_1hr/';%'/export/data/miawarac/MIAC_calibration_GROSOM-harmo/'
calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.filename];
calibrationTool.meteoFolder=['/mnt/lake/instrumentdata/meteo/zimm/Meteo/' calibrationTool.dateStr(1:4) '/'];%['/storage/lake/instrumentdata/miawarac/' dateStr(1:4) '/'];%'/home/franziska/Documents/MW/play_MIA-C_calibration/';%AB
calibrationTool.read_meteo_data =@(calibrationTool) read_meteo_data_zimmerwald(calibrationTool);
calibrationTool.filenameLookup = '/home/alistair/MIAWARA_ret/extra_files/netCDF_fields_LN2.csv';

% files for the calibration
calibrationTool.elcorr_file   = '/home/alistair/MIAWARA_ret/MIAWARA_pyarts/GROMORA-harmo/files/miawarac_elcorr.mat';
calibrationTool.antenna_file  = '/home/alistair/MIAWARA_ret/MIAWARA_pyarts/GROMORA-harmo/files/miawarac_antenna.txt';

% Defining level1a filename to read (to be adapted for other users)

calibrationTool.filenameLevel1a= [calibrationTool.level1Folder 'MIAWARA_level1a' suffix '_calib_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];


calibrationTool.checkLevel0=false;

% Log file
calibrationTool.delimiter_logfile = ';';
calibrationTool.THotUnit='degreeC';

% Function for the harmonization of the log
calibrationTool.harmonize_log=@(calibrationTool, log) harmonize_log_miawara(calibrationTool, log);  

calibrationTool.doTippingCurve = false;


calibrationTool.elevationAngleAntenna=20;
calibrationTool.elevationAngleCold=60;
calibrationTool.elevationAngleHot=180;
calibrationTool.elevationAngleRef=90;
calibrationTool.elevationAngleLN2=270;


calibrationTool.elevationAngleTolerance=15;
calibrationTool.elevationAngleHotTol = 1;
calibrationTool.elevationAngleColdTol = 1;
calibrationTool.elevationAngleRefTol = 15;
calibrationTool.elevationAngleLN2Tol = 10;

%ToCheck
calibrationTool.cycleDurationCold = 10;
calibrationTool.cycleDurationSky = 7;
calibrationTool.cycleDurationHot = 10;

calibrationTool.flipped_spectra=false;

%Other checks
calibrationTool.check_positions = false;
calibrationTool.hotTemperatureStdThreshold=0.05;

% Calibration
% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid
calibrationTool.stdAntAngleThresh = 0.5;
calibrationTool.adcOverloadThresh = 0;

% Filters for flagging "bad channels"
calibrationTool.maxStdDevTbCal = 25; 
calibrationTool.maxStdDevTbInt = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add meteo data to calibrated spectra
% calibrationTool.add_meteo_data = @(calibrationTool,correctedSpectra) add_meteo_data_unibe(calibrationTool,correctedSpectra);
% calibrationTool.add_meteo_data = @(calibrationTool, meteoData, correctedSpectra) add_meteo_data_generic(calibrationTool, meteoData, correctedSpectra);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level0 -> Level1a functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting the functions that will be used for processing this
% calibration

% Reading routine to use for the raw data
calibrationTool.read_level0=@(calibrationTool, rawFileReading) read_level0_missing(calibrationTool, rawFileReading); 

% Quality check for the raw data
calibrationTool.check_level0=@(log,rawSpectra,calibrationTool) check_level0_generic(log,rawSpectra,calibrationTool);

% Reformatting of the raw spectra into a matrix (numberOfSpectra x
% numberOfChannels)
calibrationTool.reformat_spectra=@(rawSpectra,log,calibrationTool) reformat_spectra_generic(rawSpectra,log,calibrationTool);

% Plotting some raw spectra:
calibrationTool.plot_raw_spectra=@(rawSpectra,lowerLim,upperLim,N) plot_raw_spectra_generic(rawSpectra,lowerLim,upperLim,N);


% Function to use for doing the calibration:
calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,calType) run_LN2_calibration_miawara(rawSpectra,log,calibrationTool,calType);

% Plot some calibrated spectra:
%calibrationTool.plot_calibrated_spectra=@(calibrationTool,drift,meteoData, calibratedSpectra,N) plot_spectra_generic(calibrationTool,drift,meteoData, calibratedSpectra,N);

% Function for quality check of the calibrated spectra
calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_miawara_LN2(log,calibrationTool,calibratedSpectra);

% Function saving the calibrated spectra into netCDF file

dat  = datenum(calibrationTool.dateStr);
if datenum(calibrationTool.dateStr) < datenum(2022,10,17)
    disp('date is before USRP installed')
    calibrationTool.spectrometerQuantity = 1;
    calibrationTool.numberOfSpectrometer=1;
    calibrationTool.spectrometer='AC240';
    calibrationTool.numberOfChannels = 16384;
    calibrationTool.channel_freqs = '/home/alistair/MIAWARA_ret/MIAWARA_FREQ_TST.txt';
    calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_l1a_LN2_UD(calibrationTool,log,calibratedSpectra,warningLevel0);
elseif cal_upgradeDONE
    calibrationTool.spectrometerQuantity = 3;
    calibrationTool.numberOfChannels = 49152;
    calibrationTool.numberOfSpectrometer=[1 2 3];
    calibrationTool.spectroSubChan = [16384 16384 16384];
    calibrationTool.channel_freqs = '/home/alistair/MIAWARA_ret/MIAWARA_pyarts/GROMORA-harmo/files/duel_spectrometer_freqs.dat';
    calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_l1a_LN2_UD(calibrationTool,log,calibratedSpectra,warningLevel0);
else
    calibrationTool.spectrometerQuantity = 2;
    calibrationTool.spectrometer1.spectrometer='AC240';
    calibrationTool.spectrometer1.numberOfSpectrometer=1;
    calibrationTool.spectrometer1.numberOfChannels = 16384;
    calibrationTool.spectrometer2.spectrometer='USRP';
    calibrationTool.spectrometer2.numberOfChannels = 32768;
    calibrationTool.spectrometer2.numberOfSpectrometer = 2;
    calibrationTool.save_level1a=@(calibrationTool,log,calibratedSpectra,warningLevel0) save_l1a_LN2(calibrationTool,log,calibratedSpectra,warningLevel0);
    calibrationTool.channel_freqs = '/home/alistair/MIAWARA_ret/MIAWARA_pyarts/GROMORA-harmo/files/duel_spectrometer_freqs.dat';
end
    
%Use frequencies directly from the import file
calibrationTool.read_freq_direct = true;
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


end
