function [calib, meteoData, calibrationTool] = read_level1a_daily(calibrationTool)
%==========================================================================
% NAME          | read_level1a_daily.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Function to read a level1a previously saved.
%               |
%               |
%               |
% ARGUMENTS     | INPUTS: -retrievalTool (containing the filename of the
%               | level1a file to read.
%               |
%               | OUTPUTS: - calib (part of "calibratedSpectra")
%               |          - retrievalTool
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================

% filename:
filename=calibrationTool.filenameLevel1a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate variables, directly adding the attributes
correctedSpectra.meanTime = ncread(filename,'/spectrometer1/time')';
correctedSpectra.timeUnit = ncreadatt(filename,'/spectrometer1/time','units');
correctedSpectra.timeCalendar= ncreadatt(filename,'/spectrometer1/time','calendar');

correctedSpectra.meanDateTime = datetime(correctedSpectra.meanTime+datenum(1970,1,1),'ConvertFrom','datenum');

correctedSpectra.channelID = ncread(filename,'/spectrometer1/channel_idx')';

% some variable for better identifying the time period of the measurements
correctedSpectra.year=ncread(filename,'/spectrometer1/year');
correctedSpectra.month= ncread(filename,'/spectrometer1/month');
correctedSpectra.day = ncread(filename,'/spectrometer1/day');
% ncread(filename,'/spectrometer1/timeOfDay',[calibratedSpectra.timeOfDay]);

correctedSpectra.firstSkyTime=ncread(filename,'/spectrometer1/first_sky_time')';
% ncreadatt(filename,'/spectrometer1/firstSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/spectrometer1/firstSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/spectrometer1/firstSkyTime','description','start time of the first sky measurements in this cycle');

correctedSpectra.lastSkyTime=ncread(filename,'/spectrometer1/last_sky_time')';
% ncreadatt(filename,'/spectrometer1/lastSkyTime','units',calibratedSpectra(1).meanDatetimeUnit);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','calendar',calibratedSpectra(1).calendar);
% ncreadatt(filename,'/spectrometer1/lastSkyTime','description','stop time of the first sky measurements in this cycle');

correctedSpectra.timeMin=ncread(filename,'/spectrometer1/time_min')';
correctedSpectra.tod = ncread(filename,'/spectrometer1/time_of_day')';

% % Scientific Dataset (spectrometer1,2,...)
% Calibration variable
correctedSpectra.Tb=ncread(filename,'/spectrometer1/Tb')';
correctedSpectra.stdTb=ncread(filename,'/spectrometer1/stdTb')';
correctedSpectra.freq=ncread(filename,'/spectrometer1/frequencies')';
correctedSpectra.if=ncread(filename,'/spectrometer1/intermediate_freq')';
correctedSpectra.THot=ncread(filename,'/spectrometer1/THot')';
correctedSpectra.stdTHot=ncread(filename,'/spectrometer1/stdTHot')';
correctedSpectra.TSys=ncread(filename,'/spectrometer1/TSys')';
correctedSpectra.stdTSys=ncread(filename,'/spectrometer1/stdTSys')';
correctedSpectra.meanStdTb=ncread(filename,'/spectrometer1/mean_std_Tb')';
correctedSpectra.meanAngleAntenna=ncread(filename,'/spectrometer1/mean_sky_elevation_angle')';
correctedSpectra.numHotSpectra=ncread(filename,'/spectrometer1/number_of_hot_spectra')';
correctedSpectra.numColdSpectra=ncread(filename,'/spectrometer1/number_of_cold_spectra')';
correctedSpectra.numAntSpectra=ncread(filename,'/spectrometer1/number_of_sky_spectra')';
correctedSpectra.TWindow=ncread(filename,'/spectrometer1/TWindow')';
correctedSpectra.calibrationTime=ncread(filename,'/spectrometer1/calibration_time')';
try
    correctedSpectra.noiseLevel=ncread(filename,'/spectrometer1/noise_level')';
catch ME
     warning(ME.identifier,': no noise level defined for calibration');
     correctedSpectra.noiseLevel = -9999*ones(1,length(correctedSpectra.meanTime));
end
%correctedSpectra.effectiveCalibrationTime=ncread(filename,'/spectrometer1/effectiveCalibrationTime')';

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
% Reading the flags variables
correctedSpectra.flagVector=ncread(filename,'/flags/calibration_flags')';
calibrationTool.flagVectorLength = ncreadatt(filename,'/flags/','number_of_flags');
calibrationTool.logFile.calibration_flags_meaning = [];
for i = 1:calibrationTool.flagVectorLength
    varName = ['errorCode_' num2str(i)];
    errorCodeMeaning = ncreadatt(filename,'/flags/calibration_flags',varName);
    calibrationTool.logFile.calibration_flags_meaning = [calibrationTool.logFile.calibration_flags_meaning, string(errorCodeMeaning)];
end

% calibrationTool.logFile.calibration_flags_meaning = [
%     string(ncreadatt(filename,'/flags/calibration_flags','errorCode_1')),...
%     string(ncreadatt(filename,'/flags/calibration_flags','errorCode_2')),...
%     string(ncreadatt(filename,'/flags/calibration_flags','errorCode_3')),...
%     string(ncreadatt(filename,'/flags/calibration_flags','errorCode_4')),...
%     string(ncreadatt(filename,'/flags/calibration_flags','errorCode_5')),...
%     string(ncreadatt(filename,'/flags/calibration_flags','errorCode_6'))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading the meteo variables
meteoData.dateNum = ncread(filename,'/meteo/time')';
meteoData.dateTime = datetime(meteoData.dateNum + datenum(1970,1,1),'ConvertFrom','datenum');
meteoData.air_pressure = ncread(filename,'/meteo/air_pressure')';
meteoData.air_temperature = ncread(filename,'/meteo/air_temperature')';
meteoData.relative_humidity = ncread(filename,'/meteo/relative_humidity')';
meteoData.precipitation = ncread(filename,'/meteo/precipitation')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading global attributes
calibrationTool.logFile.raw_file_warning=ncreadatt(filename,'/','raw_file_warning');
calibrationTool.logFile.comment=ncreadatt(filename,'/','comment');
calibrationTool.logFile.raw_file_comment=ncreadatt(filename,'/','raw_file_comment');
calibrationTool.logFile.rawFilename=ncreadatt(filename,'/','raw_filename');
calibrationTool.logFile.rawData=ncreadatt(filename,'/','raw_data');
calibrationTool.logFile.raw_data_software_version=ncreadatt(filename,'/','raw_data_software_version');
calibrationTool.logFile.calibration_version=ncreadatt(filename,'/','calibration_version');
calibrationTool.logFile.creation_date_level1a=ncreadatt(filename,'/','creation_date');
calibrationTool.logFile.raw_data_software_version=ncreadatt(filename,'/','raw_data_software_version');
calibrationTool.logFile.filenameLevel1a=ncreadatt(filename,'/','filename');

disp(['File read : ' filename])

for i = 1:length(correctedSpectra.meanTime)
    calib(i).Tb = correctedSpectra.Tb(i,:);
    calib(i).stdTb = correctedSpectra.stdTb(i,:);
    calib(i).dateTime = correctedSpectra.meanDateTime(i);
    calib(i).dateTime.TimeZone = calibrationTool.timeZone;
    calib(i).frequencies = correctedSpectra.freq(:)';
    calib(i).intermediate_freq = correctedSpectra.if(:)';
    calib(i).year = correctedSpectra.year(i);
    calib(i).month = correctedSpectra.month(i);
    calib(i).day = correctedSpectra.day(i);
    calib(i).THot = correctedSpectra.THot(i);
    calib(i).stdTHot = correctedSpectra.stdTHot(i);
    calib(i).TSys = correctedSpectra.TSys(i);
    calib(i).TWindow = correctedSpectra.TWindow(i);
    calib(i).meanStdTb = correctedSpectra.meanStdTb(i);
    calib(i).mean_sky_elevation_angle  =  correctedSpectra.meanAngleAntenna(i);
    calib(i).number_of_hot_spectra  =  correctedSpectra.numHotSpectra(i);
    calib(i).number_of_cold_spectra  =  correctedSpectra.numColdSpectra(i);
    calib(i).number_of_sky_spectra  =  correctedSpectra.numAntSpectra(i);
    calib(i).stdTSys = correctedSpectra.stdTSys(i);
    calib(i).calibration_time = correctedSpectra.calibrationTime(i);
    calib(i).first_sky_time = correctedSpectra.firstSkyTime(i);
    calib(i).last_sky_time = correctedSpectra.lastSkyTime(i);
    calib(i).time_min = correctedSpectra.timeMin(i);
    calib(i).time_of_day = correctedSpectra.tod(i);
    calib(i).noise_level = correctedSpectra.noiseLevel(i);
    calib(i).flags = correctedSpectra.flagVector(i,:);
end

end


