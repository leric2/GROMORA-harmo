function calib = read_level1a_daily(retrievalTool)
%==========================================================================
% NAME          | read_level1a_daily.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      |
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:
%               |
%               | OUTPUTS:
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
%
%

% filename:
filename=retrievalTool.filenameLevel1a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scientific Dataset (aquirisFFT,2,...)
% the variables linked with the calibration

correctedSpectra.Tb=ncread(filename,'/aquirisFFT/Tb')';
correctedSpectra.freq=ncread(filename,'/aquirisFFT/frequencies')';
correctedSpectra.THot=ncread(filename,'/aquirisFFT/THot')';
correctedSpectra.stdTHot=ncread(filename,'/aquirisFFT/stdTHot')';
correctedSpectra.TSys=ncread(filename,'/aquirisFFT/TSys')';
correctedSpectra.stdTSys=ncread(filename,'/aquirisFFT/stdTSys')';

correctedSpectra.calibrationTime=ncread(filename,'/aquirisFFT/calibrationTime')';
correctedSpectra.effectiveCalibrationTime=ncread(filename,'/aquirisFFT/effectiveCalibrationTime')';

% ncread(filename,'/aquirisFFT/stdTHot',[calibratedSpectra.stdTHot]);
% ncread(filename,'/aquirisFFT/TSys',[calibratedSpectra.Tsys]);
% ncread(filename,'/aquirisFFT/stdTSys',[calibratedSpectra.stdTSys]);
% ncread(filename,'/aquirisFFT/calibrationTime',[calibratedSpectra.calibrationTime]);
% ncread(filename,'/aquirisFFT/meanAngleAntenna',[calibratedSpectra.meanAngleAntenna]);
% 
% ncread(filename,'/aquirisFFT/TRoom',[calibratedSpectra.TempRoom]);
% ncread(filename,'/aquirisFFT/stdTRoom',[calibratedSpectra.stdTempRoom]);
% ncread(filename,'/aquirisFFT/TOut',[calibratedSpectra.TempOut]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the flags variables
% ncread(filename,'/flags/time',[calibratedSpectra.meanDatetime]);
% ncreadatt(filename,'/flags/time','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/flags/time','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/flags/time','description','mean time of the measurements for this cycle');
% 
% ncread(filename,'/flags/flags',1:length(calibratedSpectra(1).errorVector));
% ncread(filename,'/flags/calibration_flags',errorCalib');
% 

% Writing the aquirisFFT variable
% Coordinate variables, directly adding the attributes
correctedSpectra.meanTime = ncread(filename,'/aquirisFFT/time')';
correctedSpectra.timeUnit = ncreadatt(filename,'/aquirisFFT/time','units');
correctedSpectra.timeCalendar= ncreadatt(filename,'/aquirisFFT/time','calendar');

correctedSpectra.channelID = ncread(filename,'/aquirisFFT/channel_idx')';

% some variable for better identifying the time period of the measurements
correctedSpectra.year=ncread(filename,'/aquirisFFT/year');
correctedSpectra.month= ncread(filename,'/aquirisFFT/month');
correctedSpectra.day = ncread(filename,'/aquirisFFT/day');
% ncread(filename,'/aquirisFFT/timeOfDay',[calibratedSpectra.timeOfDay]);

correctedSpectra.firstSkyTime=ncread(filename,'/aquirisFFT/firstSkyTime')';  
% ncreadatt(filename,'/aquirisFFT/firstSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/aquirisFFT/firstSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/aquirisFFT/firstSkyTime','description','start time of the first sky measurements in this cycle');

correctedSpectra.lastSkyTime=ncread(filename,'/aquirisFFT/lastSkyTime')';
% ncreadatt(filename,'/aquirisFFT/lastSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/aquirisFFT/lastSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/aquirisFFT/lastSkyTime','description','stop time of the first sky measurements in this cycle');

disp(['File read : ' filename])

for i = 1:length(correctedSpectra.meanTime)
    calib(i).Tb=correctedSpectra.Tb(i,:);
    calib(i).meanDatetime=correctedSpectra.meanTime(i);
    calib(i).freq=correctedSpectra.freq(:)';
    calib(i).year=correctedSpectra.year(i);
    calib(i).month=correctedSpectra.month(i);
    calib(i).day=correctedSpectra.day(i);
    calib(i).THot=correctedSpectra.THot(i);
    calib(i).stdTHot=correctedSpectra.stdTHot(i);
    calib(i).TSys=correctedSpectra.TSys(i);
    calib(i).stdTSys=correctedSpectra.stdTSys(i);
    calib(i).calibrationTime=correctedSpectra.calibrationTime(i);
    calib(i).firstSkyTime=correctedSpectra.firstSkyTime(i);
    calib(i).lastSkyTime=correctedSpectra.lastSkyTime(i);
    
end
