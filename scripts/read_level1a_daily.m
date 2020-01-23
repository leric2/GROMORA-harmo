function correctedSpectra = read_level1a_daily(retrievalTool)
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
% Scientific Dataset (spectrometer1,2,...)
% the variables linked with the calibration

correctedSpectra.Tb=ncread(filename,'/spectrometer1/Tb')';
correctedSpectra.THot=ncread(filename,'/spectrometer1/THot')';
correctedSpectra.stdTHot=ncread(filename,'/spectrometer1/stdTHot')';
correctedSpectra.TSys=ncread(filename,'/spectrometer1/TSys')';
correctedSpectra.stdTSys=ncread(filename,'/spectrometer1/stdTSys')';

correctedSpectra.calibrationTime=ncread(filename,'/spectrometer1/calibrationTime')';
correctedSpectra.effectiveCalibrationTime=ncread(filename,'/spectrometer1/effectiveCalibrationTime')';

% ncread(filename,'/spectrometer1/stdTHot',[calibratedSpectra.stdTHot]);
% ncread(filename,'/spectrometer1/TSys',[calibratedSpectra.Tsys]);
% ncread(filename,'/spectrometer1/stdTSys',[calibratedSpectra.stdTSys]);
% ncread(filename,'/spectrometer1/calibrationTime',[calibratedSpectra.calibrationTime]);
% ncread(filename,'/spectrometer1/meanAngleAntenna',[calibratedSpectra.meanAngleAntenna]);
% 
% ncread(filename,'/spectrometer1/TRoom',[calibratedSpectra.TempRoom]);
% ncread(filename,'/spectrometer1/stdTRoom',[calibratedSpectra.stdTempRoom]);
% ncread(filename,'/spectrometer1/TOut',[calibratedSpectra.TempOut]);

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

% Writing the spectrometer1 variable
% Coordinate variables, directly adding the attributes
correctedSpectra.meanTime = ncread(filename,'/spectrometer1/time')';
correctedSpectra.timeUnit = ncreadatt(filename,'/spectrometer1/time','units');
correctedSpectra.timeCalendar= ncreadatt(filename,'/spectrometer1/time','calendar');

correctedSpectra.channelID = ncread(filename,'/spectrometer1/channel_idx')';

% some variable for better identifying the time period of the measurements
% ncread(filename,'/spectrometer1/year',int64([calibratedSpectra.year]));
% ncread(filename,'/spectrometer1/month',int64([calibratedSpectra.month]));
% ncread(filename,'/spectrometer1/day',int64([calibratedSpectra.day]));
% ncread(filename,'/spectrometer1/timeOfDay',[calibratedSpectra.timeOfDay]);

correctedSpectra.firstSkyTime=ncread(filename,'/spectrometer1/firstSkyTime')';  
% ncreadatt(filename,'/spectrometer1/firstSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/spectrometer1/firstSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/spectrometer1/firstSkyTime','description','start time of the first sky measurements in this cycle');

correctedSpectra.lastSkyTime=ncread(filename,'/spectrometer1/lastSkyTime')';
% ncreadatt(filename,'/spectrometer1/lastSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','description','stop time of the first sky measurements in this cycle');


disp(['File read : ' filename])
end
