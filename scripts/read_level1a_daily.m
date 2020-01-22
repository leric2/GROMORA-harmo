function read_level1a_daily(retrievalTool)
%==========================================================================
% NAME          | read.m
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
% Saving level1a into DAILY netCDF file
% We are using the netcdf package because it offers far more flexibility
% for the writing.
locationLevel1a=retrievalTool.level1Folder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here LOOOP in all calibration cycle and write every variable into daily vector
% of matrices

%initialize some variable (only the matrices)
channelId=int64(ones(length(calibratedSpectra),retrievalTool.numberOfChannels)*NaN);
Tb=ones(length(calibratedSpectra),retrievalTool.numberOfChannels)*NaN;

error=int64(ones(length(calibratedSpectra),length(calibratedSpectra(1).errorVector))*NaN);
errorCalib=int64(ones(length(calibratedSpectra),length(calibratedSpectra(1).errorVector))*NaN);

for t = 1:length(calibratedSpectra)
    channelId(t,:)=1:retrievalTool.numberOfChannels;
    Tb(t,:)=calibratedSpectra(t).Tb;
    
    error(t,:)=1:length(calibratedSpectra(1).errorVector);
    errorCalib(t,:)=calibratedSpectra(t).errorVector;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now write daily level 1a file
filename=retrievalTool.filenameLevel1a;
%filename=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(t) '.nc'];
%title=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(i)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Scientific Dataset (spectrometer1,2,...)
%%%%%%%%%%%%%%%%%
% the variables linked with the calibration
ncread(filename,'/spectrometer1/effectiveCalibrationTime',[calibratedSpectra.effectiveCalibrationTime]);

calibratedSpectra.Tb=ncread(filename,'/spectrometer1/Tb')';
calibratedSpectra.THot=ncread(filename,'/spectrometer1/THot')';
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
ncread(filename,'/spectrometer1/time',[calibratedSpectra.meanDatetime]);
ncreadatt(filename,'/spectrometer1/time','units',calibratedSpectra(1).meanDatetimeUnit);
ncreadatt(filename,'/spectrometer1/time','calendar',calibratedSpectra(1).calendar);
ncreadatt(filename,'/spectrometer1/time','description','mean time of the measurements for this cycle');

ncread(filename,'/spectrometer1/channel_idx',1:retrievalTool.numberOfChannels);

% some variable for better identifying the time period of the measurements
% ncread(filename,'/spectrometer1/year',int64([calibratedSpectra.year]));
% ncread(filename,'/spectrometer1/month',int64([calibratedSpectra.month]));
% ncread(filename,'/spectrometer1/day',int64([calibratedSpectra.day]));
% ncread(filename,'/spectrometer1/timeOfDay',[calibratedSpectra.timeOfDay]);

ncread(filename,'/spectrometer1/firstSkyTime',[calibratedSpectra.datetimeStart]);  
% ncreadatt(filename,'/spectrometer1/firstSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/spectrometer1/firstSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/spectrometer1/firstSkyTime','description','start time of the first sky measurements in this cycle');

ncread(filename,'/spectrometer1/lastSkyTime',[calibratedSpectra.datetimeStop]);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','description','stop time of the first sky measurements in this cycle');


disp(['File read : ' filename])
end
