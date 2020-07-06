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

% Name of the instrument
calibrationTool.instrumentName=instrumentName;
calibrationTool.dateStr=dateStr;

% Valid properties for all instruments
calibrationTool.lightSpeed=299792458; % [m/s] 
calibrationTool.bytesPerValue=4;
calibrationTool.binaryType='ieee-be';

   
calibrationTool.h=6.62606876e-34; % [J/s]
calibrationTool.kb=1.38065e-23;    % [J/K]

%calibrationTool.hotTemperatureStdThreshold=0.05;

% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid
%calibrationTool.minNumberOfIndicePerCycle=5;

%calibrationTool.numberOfAquisitionSpectraHot=60000;
%calibrationTool.numberOfAquisitionSpectraAntenna=120000;
%calibrationTool.numberOfAquisitionSpectraCold=120000;

%calibrationTool.TSysThresh=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing default parameters for each instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROMOS
switch instrumentName
    case 'GROMOS'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Meta data
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:
        calibrationTool.numberOfSpectrometer=1;
        calibrationTool.spectrometer='AC240';
        
        calibrationTool.samplingRateFFTS=2000;
        
        calibrationTool.numberOfChannels=32768;
        
        % Local oscillators information
        calibrationTool.fLO1=1.45875e11;
        calibrationTool.fLO2=3.6e9;
        
        % This one should correspond to the DC channel
        calibrationTool.LOFreqTot=calibrationTool.fLO1-calibrationTool.fLO2-0.5e9;
        
        %calibrationTool.DCChannel=16384;
        calibrationTool.instrumentBandwidth=1e9;
        %calibrationTool.LOFreq=1.45875e11;
  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for flagging
        calibrationTool.TSysCenterTh=2750;
        calibrationTool.TSysThresh=100;
        calibrationTool.stdTSysThresh=8;
        
        calibrationTool.THotTh=313.9;
        calibrationTool.THotAbsThresh=2;
        
        
        calibrationTool.numberOfTippingCurveExpected=48;
        calibrationTool.toleranceTippingCurves=2;
        calibrationTool.elevationAngleAntenna=40;
        calibrationTool.elevationAngleCold=-84;
        % -84 before 2020 ?????
        %calibrationTool.elevationAngleCold=-84;
        calibrationTool.elevationAngleHot=160;
        
        calibrationTool.elevationAngleTolerance=5;
        % Considering the expected number of tipping curve:
        calibrationTool.numberOfCyclesExpected=1500;
        calibrationTool.toleranceNumberCycles=0.01*calibrationTool.numberOfCyclesExpected;
        calibrationTool.tippingSize=27;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Corrections
        calibrationTool.tWindow=0.99;
        
        calibrationTool.checkLevel0=true;

        calibrationTool.flipped_spectra=true;
        calibrationTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Log file from the instrument:
        calibrationTool.THotUnit='degreeC';
        % Function for the harmonization of the log
        calibrationTool.harmonize_log=@(log) harmonize_log_gromos(log);  
        
        calibrationTool.indiceCold=0;
        calibrationTool.indiceAntenna=1;
        calibrationTool.indiceHot=2;
        
        % Calibration outlier management
        calibrationTool.threshNumRawSpectraHot=0.05*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraCold=0.05*calibrationTool.numberOfChannels;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOMORA    
    case 'SOMORA'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Meta data
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:
        calibrationTool.numberOfSpectrometer=1;
        calibrationTool.spectrometer='AC240';
        
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for flagging
        
        calibrationTool.TSysCenterTh=2750;
        calibrationTool.TSysThresh=100;
        calibrationTool.stdTSysThresh=8;
        
        calibrationTool.THotTh=311.1;
        calibrationTool.THotAbsThresh=2;
        
        calibrationTool.checkLevel0=true;
        
        calibrationTool.numberOfTippingCurveExpected=4;
        calibrationTool.toleranceTippingCurves=2;
        
        
        calibrationTool.numberOfCyclesExpected=3940;
        calibrationTool.toleranceNumberCycles=0.01*calibrationTool.numberOfCyclesExpected;
        
        calibrationTool.tippingSize=5;
        calibrationTool.elevationAngleAntenna=38;
        calibrationTool.elevationAngleCold=-90;
        calibrationTool.elevationAngleHot=180;
        calibrationTool.elevationAngleTolerance=5;

        calibrationTool.flipped_spectra=false;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Log file from the instrument:
        % The log needs to be harmonized, for now taking gromos as basis
        calibrationTool.harmonize_log=@(log) harmonize_log_somora(log);
        calibrationTool.THotUnit='K';
        
        calibrationTool.indiceCold=0;
        calibrationTool.indiceAntenna=1;
        calibrationTool.indiceHot=2;
        
        % Calibration outlier management
        calibrationTool.threshNumRawSpectraHot=0.05*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraCold=0.05*calibrationTool.numberOfChannels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOPI       
    case 'mopi5'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Meta data
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:

        calibrationTool.numberOfChannels=16384;

        calibrationTool.instrumentBandwidth=1e9;
        calibrationTool.LOFreq=1.45875e11;

        calibrationTool.indiceCold=2;
        calibrationTool.indiceAntenna=5;
        calibrationTool.indiceHot=1;
        
        calibrationTool.checkLevel0=false;
        
        %calibrationTool.numberOfTippingCurveExpected=48;
        %calibrationTool.toleranceTippingCurves=2;
        calibrationTool.elevationAngleAntenna=140;
        calibrationTool.elevationAngleCold=265;
        calibrationTool.elevationAngleHot=85;
        
        calibrationTool.elevationAngleTolerance=5;
        % Considering the expected number of tipping curve:
        %calibrationTool.numberOfCyclesExpected=1500;
        %calibrationTool.toleranceNumberCycles=0.01*calibrationTool.numberOfCyclesExpected;
        calibrationTool.tippingSize=27;
        calibrationTool.flipped_spectra=false;
        calibrationTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);
        %calibrationTool.THotUnit='degreeC';
        % Harmonizing mopi5 log (for units)
        calibrationTool.harmonize_log=@(log) harmonize_log_mopi5(log);  
        
        % Calibration outlier management
        calibrationTool.threshNumRawSpectraHot=0.05*calibrationTool.numberOfChannels;
        calibrationTool.threshNumRawSpectraCold=0.05*calibrationTool.numberOfChannels;
        
        calibrationTool.rawFileFolder=['/mnt/instrumentdata/mopi5/' dateStr(1:4) '/'];
        %calibrationTool.rawFileFolder=['/scratch/mopi_rawData/'];
        calibrationTool.level1Folder='/home/esauvageat/Documents/MOPI5/Level1/';
        
        calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.instrumentName,'_', calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        
        calibrationTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        calibrationTool.observationFreq=110;
        
        calibrationTool.calibrationTime=60;
        
        %calibrationTool.fLO1=1.49275e11;
        %calibrationTool.fLO2=5.6e9;
        %calibrationTool.fLO3=2e9;
        
        % This one should correspond to the DC channel
        calibrationTool.LOFreqTot=1.10e11;
        calibrationTool.DCChannel=1; %=Nchannel/2 ??
        

        calibrationTool.ffts_model=1;
        calibrationTool.numberOfSpectrometer=4;
        S  = {'USRP-A', 'USRP-B','U5303', 'AC240'};
        calibrationTool.spectrometer=S{calibrationTool.ffts_model};
        
        FS = [200 20  3200 2000]; % sampling rates in MHz 
        calibrationTool.FS=FS(calibrationTool.ffts_model);
        calibrationTool.read_level0=@(calibrationTool) mopi5_read(calibrationTool); 
        
        calibrationTool.filenameLevel1a=['/home/esauvageat/Documents/MOPI5/Level1/mopi5_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];
        
        calibrationTool.calibrate=@(rawSpectra,log,calibrationTool,TCold,calType) calibrate_mopi5(rawSpectra,log,calibrationTool,TCold,calType);
        calibrationTool.check_calibrated=@(log,calibrationTool,calibratedSpectra) check_calibrated_mopi5(log,calibrationTool,calibratedSpectra);

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
        calibrationTool.harmonize_log=@(log) harmonize_log_miawarac(log);  
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
        calibrationTool.rawFileFolder='/home/franziska/Documents/MW/play_MIA-C_calibration/';%['/mnt/instrumentdata/miawarac/' dateStr(1:4) '/'];
        %calibrationTool.rawFileFolder=['/scratch/mopi_rawData/'];
        calibrationTool.level1Folder='/home/franziska/Documents/MW/play_MIA-C_calibration/';
        
        calibrationTool.file=[calibrationTool.rawFileFolder,calibrationTool.instrumentName,'_', calibrationTool.dateStr(1:4) '_' calibrationTool.dateStr(6:7) '_' calibrationTool.dateStr(9:10)];
        
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

        calibrationTool.read_level0=@(calibrationTool) read_level0_missing(calibrationTool); 
        
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
        % ...
    case 'SOMORA'
        % window transmission 
        if timeNumber<datenum(2007,07,03) %attn anciennement  733957 cad 3 juillet 2009 mais FAUX
            calibrationTool.tWindow = 0.987; %value from 2002 to 3/7/2007 (att. window has been changed 8/2/2006 but new t not measured) nb: ema 9/11/2007
        elseif (timeNumber>= datenum(2007,07,03) && timeNumber<datenum(2010,6,26)) %ie 734315
            calibrationTool.tWindow = 0.979; %value from 3/7/2007 to 25/6/2010   nb: ema 9/11/2007; calc sur la base de R-J
        elseif (timeNumber>= datenum(2010,06,25) && timeNumber<datenum(2012,03,02)) %ie 734931 2/3/2012
            calibrationTool.tWindow = 0.985; %value from 26/6/2010 to 2/3/2012 %attn anciennement 0.976 ???
        elseif (timeNumber>= datenum(2012,03,02) && timeNumber<datenum(2014,06,16)) %ie 16/06/2014
            calibrationTool.tWindowt = 0.999; %value from 3/3/2012 to 15/06/2014
        elseif (timeNumber>= datenum(2014,06,16) && timeNumber<datenum(2015,09,08))
            calibrationTool.tWindow = 1.0105; %value from 16/6/2014 to 07/09/2015
        elseif (timeNumber>= datenum(2015,09,09) && timeNumber<datenum(2017,05,30))
            calibrationTool.tWindow = 1.009; %value from 8/09/2015 to 29/05/2017
        elseif (timeNumber>= datenum(2017,05,31) && timeNumber<datenum(2018,08,21))
            calibrationTool.tWindow = 0.9922; %value from 30/05/2017 to 20/08/2018
        else
            calibrationTool.tWindow = 0.9980; %value since 21/08/2018
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
