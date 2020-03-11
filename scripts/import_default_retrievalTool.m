function retrievalTool = import_default_retrievalTool(instrumentName,dateStr)
%==========================================================================
% NAME          | import_default_retrievalTool(instrumentName)
% TYPE          | Function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: * instrumentName: name of the instrument as a
%               | string (ex: 'GROMOS')
%               |
%               | OUTPUTS: * retrievalTool: the default toolbox for
%               | launching a retrieval for this instrument.
%               | 
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
retrievalTool = struct();

% Name of the instrument
retrievalTool.instrumentName=instrumentName;
retrievalTool.dateStr=dateStr;

% Valid properties for all instruments
retrievalTool.lightSpeed=299792458; % [m/s] 
retrievalTool.bytesPerValue=4;
retrievalTool.binaryType='ieee-be';

   
retrievalTool.h=6.62606876e-34; % [J/s]
retrievalTool.kb=1.38065e-23;    % [J/K]

retrievalTool.hotTemperatureStdThreshold=0.05;

% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid
retrievalTool.minNumberOfIndicePerCycle=5;

retrievalTool.numberOfAquisitionSpectraHot=60000;
retrievalTool.numberOfAquisitionSpectraAntenna=120000;
retrievalTool.numberOfAquisitionSpectraCold=120000;

retrievalTool.TSysThresh=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing default parameters for each instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROMOS
switch string(instrumentName)
    case 'GROMOS'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Meta data
        retrievalTool.dataLocation='BERN';
        retrievalTool.PI_NAME='Murk;Axel';
        retrievalTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
        retrievalTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
        retrievalTool.PI_EMAIL='axel.murk@iap.unibe.ch'; 
        retrievalTool.dataSource='MWR.O3_UBERN';
        
        % Geolocation
        retrievalTool.lon=7.44;
        retrievalTool.lat=46.95;
        retrievalTool.altitude=560;
        retrievalTool.azimuthAngle=45;
        
        % Observation frequency
        retrievalTool.observationFreq=1.4217504e11;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:
        retrievalTool.numberOfSpectrometer=1;
        retrievalTool.spectrometer='AC240';
        
        retrievalTool.samplingRateFFTS=2000;
        
        retrievalTool.numberOfChannels=32768;
        
        % Local oscillators information
        retrievalTool.fLO1=1.45875e11;
        retrievalTool.fLO2=3.6e9;
        
        % This one should correspond to the DC channel
        retrievalTool.LOFreqTot=retrievalTool.fLO1-retrievalTool.fLO2-0.5e9;
        
        %retrievalTool.DCChannel=16384;
        retrievalTool.instrumentBandwidth=1e9;
        %retrievalTool.LOFreq=1.45875e11;
  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for flagging
        retrievalTool.TSysCenterTh=2750;
        retrievalTool.TSysThresh=100;
        retrievalTool.stdTSysThresh=8;
        
        retrievalTool.THotTh=313.9;
        retrievalTool.THotAbsThresh=2;
        
        
        retrievalTool.numberOfTippingCurveExpected=48;
        retrievalTool.toleranceTippingCurves=2;
        retrievalTool.elevationAngleAntenna=40;
        retrievalTool.elevationAngleCold=-84;
        retrievalTool.elevationAngleHot=160;
        
        retrievalTool.elevationAngleTolerance=5;
        % Considering the expected number of tipping curve:
        retrievalTool.numberOfCyclesExpected=1500;
        retrievalTool.toleranceNumberCycles=0.01*retrievalTool.numberOfCyclesExpected;
        retrievalTool.tippingSize=27;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Corrections
        retrievalTool.tWindow=0.99;
        
        retrievalTool.checkLevel0=true;

        retrievalTool.flipped_spectra=true;
        retrievalTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Log file from the instrument:
        retrievalTool.THotUnit='degreeC';
        % Function for the harmonization of the log
        retrievalTool.harmonize_log=@(log) harmonize_log_gromos(log);  
        
        retrievalTool.indiceCold=0;
        retrievalTool.indiceAntenna=1;
        retrievalTool.indiceHot=2;
        
        % Calibration outlier management
        retrievalTool.threshNumRawSpectraHot=0.05*retrievalTool.numberOfChannels;
        retrievalTool.threshNumRawSpectraCold=0.05*retrievalTool.numberOfChannels;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOMORA    
    case 'SOMORA'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Meta data
        retrievalTool.dataLocation='PAYERNE';
        retrievalTool.PI_NAME='';
        retrievalTool.PI_AFFILIATION='Swiss Meteorological Institute;MCH';
        retrievalTool.PI_ADDRESS='';
        retrievalTool.PI_EMAIL=''; 
        retrievalTool.dataSource='MWR.O3_MCH';
        
        % Geolocation
        retrievalTool.lon=6.95;
        retrievalTool.lat=46.82;
        retrievalTool.altitude=491;
        retrievalTool.azimuthAngle=34;
        
        % Observation frequency
        retrievalTool.observationFreq=1.4217504e11;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:
        retrievalTool.numberOfSpectrometer=1;
        retrievalTool.spectrometer='AC240';
        
        retrievalTool.fLO1=1.49275e11;
        retrievalTool.fLO2=5.6e9;
        retrievalTool.fLO3=2e9;
        
        % This one should correspond to the DC channel
        retrievalTool.LOFreqTot=retrievalTool.fLO1-retrievalTool.fLO2-retrievalTool.fLO3;
        %retrievalTool.DCChannel=1; %=Nchannel/2 ??

        retrievalTool.numberOfChannels=16384;

        retrievalTool.instrumentBandwidth=1e9;
        retrievalTool.LOFreq=1.4927504e11; %TOCHECK
        retrievalTool.samplingRateFFTS=2000;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for flagging
        
        retrievalTool.TSysCenterTh=2750;
        retrievalTool.TSysThresh=100;
        retrievalTool.stdTSysThresh=8;
        
        retrievalTool.THotTh=311.1;
        retrievalTool.THotAbsThresh=2;
        
        retrievalTool.checkLevel0=true;
        
        retrievalTool.numberOfTippingCurveExpected=4;
        retrievalTool.toleranceTippingCurves=2;
        
        
        retrievalTool.numberOfCyclesExpected=3940;
        retrievalTool.toleranceNumberCycles=0.01*retrievalTool.numberOfCyclesExpected;
        
        retrievalTool.tippingSize=5;
        retrievalTool.elevationAngleAntenna=38;
        retrievalTool.elevationAngleCold=-90;
        retrievalTool.elevationAngleHot=180;
        retrievalTool.elevationAngleTolerance=5;

        retrievalTool.flipped_spectra=false;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Log file from the instrument:
        % The log needs to be harmonized, for now taking gromos as basis
        retrievalTool.harmonize_log=@(log) harmonize_log_somora(log);
        retrievalTool.THotUnit='K';
        
        retrievalTool.indiceCold=0;
        retrievalTool.indiceAntenna=1;
        retrievalTool.indiceHot=2;
        
        % Calibration outlier management
        retrievalTool.threshNumRawSpectraHot=0.05*retrievalTool.numberOfChannels;
        retrievalTool.threshNumRawSpectraCold=0.05*retrievalTool.numberOfChannels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOPI       
    case 'mopi5'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Meta data
        retrievalTool.dataLocation='BERN';
        retrievalTool.PI_NAME='Murk;Axel';
        retrievalTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
        retrievalTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
        retrievalTool.PI_EMAIL='axel.murk@iap.unibe.ch'; 
        retrievalTool.dataSource='';
        
        % Geolocation
        retrievalTool.lon=7.44;
        retrievalTool.lat=46.95;
        retrievalTool.altitude=560;
        retrievalTool.azimuthAngle=45;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spectrometer data:

        retrievalTool.numberOfChannels=16384;

        retrievalTool.instrumentBandwidth=1e9;
        retrievalTool.LOFreq=1.45875e11;

        retrievalTool.indiceCold=2;
        retrievalTool.indiceAntenna=5;
        retrievalTool.indiceHot=1;
        
        retrievalTool.checkLevel0=false;
        
        %retrievalTool.numberOfTippingCurveExpected=48;
        %retrievalTool.toleranceTippingCurves=2;
        retrievalTool.elevationAngleAntenna=140;
        retrievalTool.elevationAngleCold=265;
        retrievalTool.elevationAngleHot=85;
        
        retrievalTool.elevationAngleTolerance=5;
        % Considering the expected number of tipping curve:
        %retrievalTool.numberOfCyclesExpected=1500;
        %retrievalTool.toleranceNumberCycles=0.01*retrievalTool.numberOfCyclesExpected;
        retrievalTool.tippingSize=27;
        retrievalTool.flipped_spectra=false;
        retrievalTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);
        %retrievalTool.THotUnit='degreeC';
        % Harmonizing mopi5 log (for units)
        retrievalTool.harmonize_log=@(log) harmonize_log_mopi5(log);  
        
        % Calibration outlier management
        retrievalTool.threshNumRawSpectraHot=0.05*retrievalTool.numberOfChannels;
        retrievalTool.threshNumRawSpectraCold=0.05*retrievalTool.numberOfChannels;
        
        retrievalTool.rawFileFolder=['/mnt/instrumentdata/mopi5/' dateStr(1:4) '/'];
        %retrievalTool.rawFileFolder=['/scratch/mopi_rawData/'];
        retrievalTool.level1Folder='/home/esauvageat/Documents/MOPI5/Level1/';
        
        retrievalTool.file=[retrievalTool.rawFileFolder,retrievalTool.instrumentName,'_', retrievalTool.dateStr(1:4) retrievalTool.dateStr(6:7) retrievalTool.dateStr(9:10)];
        
        retrievalTool.meteoFolder='/mnt/instrumentdata/meteo/exwi/meteo/';
        retrievalTool.observationFreq=110;
        
        retrievalTool.calibrationTime=60;
        
        %retrievalTool.fLO1=1.49275e11;
        %retrievalTool.fLO2=5.6e9;
        %retrievalTool.fLO3=2e9;
        
        % This one should correspond to the DC channel
        retrievalTool.LOFreqTot=1.10e11;
        retrievalTool.DCChannel=1; %=Nchannel/2 ??
        

        retrievalTool.ffts_model=1;
        retrievalTool.numberOfSpectrometer=4;
        S  = {'USRP-A', 'USRP-B','U5303', 'AC240'};
        retrievalTool.spectrometer=S{retrievalTool.ffts_model};
        
        FS = [200 20  3200 2000]; % sampling rates in MHz 
        retrievalTool.FS=FS(retrievalTool.ffts_model);
        retrievalTool.read_level0=@(retrievalTool) mopi5_read(retrievalTool); 
        
        retrievalTool.filenameLevel1a=['/home/esauvageat/Documents/MOPI5/Level1/mopi5_level1a_' retrievalTool.spectrometer '_' retrievalTool.dateStr '.nc'];
        
        retrievalTool.calibrate=@(rawSpectra,log,retrievalTool,TCold,calType) calibrate_mopi5(rawSpectra,log,retrievalTool,TCold,calType);
        retrievalTool.check_calibrated=@(log,retrievalTool,calibratedSpectra) check_calibrated_mopi5(log,retrievalTool,calibratedSpectra);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter varying with time for the instruments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conversion of datestr to datenum:
timeNumber=datenum(str2num(dateStr(1:4)),str2num(dateStr(6:7)),str2num(dateStr(9:10)));

switch string(instrumentName)
    case 'GROMOS'
        % ...
    case 'SOMORA'
        % window transmission 
        if timeNumber<datenum(2007,07,03) %attn anciennement  733957 cad 3 juillet 2009 mais FAUX
            retrievalTool.tWindow = 0.987; %value from 2002 to 3/7/2007 (att. window has been changed 8/2/2006 but new t not measured) nb: ema 9/11/2007
        elseif (timeNumber>= datenum(2007,07,03) && timeNumber<datenum(2010,6,26)) %ie 734315
            retrievalTool.tWindow = 0.979; %value from 3/7/2007 to 25/6/2010   nb: ema 9/11/2007; calc sur la base de R-J
        elseif (timeNumber>= datenum(2010,06,25) && timeNumber<datenum(2012,03,02)) %ie 734931 2/3/2012
            retrievalTool.tWindow = 0.985; %value from 26/6/2010 to 2/3/2012 %attn anciennement 0.976 ???
        elseif (timeNumber>= datenum(2012,03,02) && timeNumber<datenum(2014,06,16)) %ie 16/06/2014
            retrievalTool.tWindowt = 0.999; %value from 3/3/2012 to 15/06/2014
        elseif (timeNumber>= datenum(2014,06,16) && timeNumber<datenum(2015,09,08))
            retrievalTool.tWindow = 1.0105; %value from 16/6/2014 to 07/09/2015
        elseif (timeNumber>= datenum(2015,09,09) && timeNumber<datenum(2017,05,30))
            retrievalTool.tWindow = 1.009; %value from 8/09/2015 to 29/05/2017
        elseif (timeNumber>= datenum(2017,05,31) && timeNumber<datenum(2018,08,21))
            retrievalTool.tWindow = 0.9922; %value from 30/05/2017 to 20/08/2018
        else
            retrievalTool.tWindow = 0.9980; %value since 21/08/2018
        end
end
