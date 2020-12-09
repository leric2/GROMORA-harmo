function calibrationTool = save_level1a_daily(calibrationTool,logFile,calibratedSpectra,warningLevel0)
%==========================================================================
% NAME          | save_level1a_daily.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Function saving the calibrated spectra for a complete day
%               | in the form of netCDF4 file. The file contains one entry
%               | per calibration cylce (typically 10 minutes) on a time 
%               | coordinate as well as all necessery metadata to identify
%               | and re-use the data.
%               |
%               |
% ARGUMENTS     | INPUTS: - calibrationTool
%               |         - logFile
%               |         - calibratedSpectra
%               |         - warningLevel0
%               |         
%               | OUTPUTS: - saving netCDF4 level1a file
%               |          - calibrationTool with add field
%               |
% SAVE          | Level1a netCDF4 containing the following groups:
%               |   - '/' for global attributes
%               |   - '/spectrometer1/' for integration
%               |   - '/flags/' for
%==========================================================================

% Filename and location for DAILY netCDF file
if calibrationTool.calType=="debug"
    filename=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '_debug.nc'];
else
    filename=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];
end
calibrationTool.filenameLevel1a=filename;

% Rewrite if already existing
if isfile(filename)
    delete(filename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scientific Dataset (spectrometer1,2,...)

% Coordinates variable (enable 'netcdf4' format)
nccreate(filename,'/spectrometer1/time','Dimensions',{'time',Inf},'Datatype','double','Format','netcdf4');
nccreate(filename,'/spectrometer1/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','int64','FillValue',-9999)

%%%%%%%%%%%%%%%%%
% Some variables to help geolocate the file:
% these are scalar variables as they do not vary in time
nccreate(filename,'/spectrometer1/lat','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/lon','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/alt','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/azimuth_angle','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);

% time variables
nccreate(filename,'/spectrometer1/MJD2K','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/year','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/month','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/day','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/time_of_day','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/first_sky_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/last_sky_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/time_min','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);

%%%%%%%%%%%%%%%%%
% Calibration variables   
%nccreate(filename,'/spectrometer1/effectiveCalibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/Tb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/stdTb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/frequencies','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/intermediate_freq','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/mean_std_Tb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_sky_elevation_angle','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/TRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TWindow','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/noise_level','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_of_hot_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_cold_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_sky_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)

nccreate(filename,'/spectrometer1/mean_hot_counts','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags dataset

% Coordinates time and flags
nccreate(filename,'/flags/time','Dimensions',{'time',Inf},'Datatype','double')

% if the error vector does not exist, we replace it with a scalar NaN
if isfield(calibratedSpectra,'errorVector')
    lenErrorVect=length(calibratedSpectra(1).errorVector);
else
    lenErrorVect=1;
end
nccreate(filename,'/flags/flags','Dimensions',{'flags',lenErrorVect},'Datatype','int64')

% We input a (1xerrorVectorSize) int vector to identify the errors
nccreate(filename,'/flags/calibration_flags','Dimensions',{'flags',lenErrorVect,'time',Inf},'Datatype','int64','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
nccreate(filename,'/meteo/time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/air_pressure','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/air_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/relative_humidity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/precipitation','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scientific Dataset (spectrometer1,2,...)

% Coordinate variables, directly with their attributes
ncwrite(filename,'/spectrometer1/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/spectrometer1/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/time','description','mean time of THE BEGINNING of all antenna measurements for this cycle');

if isfield(calibratedSpectra, 'meanDatetimeMJD2K')
    ncwrite(filename,'/spectrometer1/MJD2K',[calibratedSpectra.meanDatetimeMJD2K]);
    ncwriteatt(filename,'/spectrometer1/MJD2K','units','MJD2K');
    ncwriteatt(filename,'/spectrometer1/MJD2K','calendar','Julian');
    ncwriteatt(filename,'/spectrometer1/MJD2K','description','mean time of THE BEGINNING of all antenna measurements for this cycle');
end
ncwrite(filename,'/spectrometer1/channel_idx',1:calibrationTool.numberOfChannels);

%%%%%%%%%%%%%%%%%
% Geolocation variables
ncwrite(filename,'/spectrometer1/lat',ones(length(calibratedSpectra),1)*calibrationTool.lat);
ncwrite(filename,'/spectrometer1/lon',ones(length(calibratedSpectra),1)*calibrationTool.lon);
ncwrite(filename,'/spectrometer1/alt',ones(length(calibratedSpectra),1)*calibrationTool.altitude);
ncwrite(filename,'/spectrometer1/azimuth_angle',ones(length(calibratedSpectra),1)*calibrationTool.azimuthAngle);

% Time variables
ncwrite(filename,'/spectrometer1/year',int64([calibratedSpectra.year]));
ncwrite(filename,'/spectrometer1/month',int64([calibratedSpectra.month]));
ncwrite(filename,'/spectrometer1/day',int64([calibratedSpectra.day]));
ncwrite(filename,'/spectrometer1/time_of_day',[calibratedSpectra.timeOfDay]);

ncwrite(filename,'/spectrometer1/first_sky_time',[calibratedSpectra.firstSkyTime]);  
ncwriteatt(filename,'/spectrometer1/first_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/first_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/first_sky_time','description','start time of the first sky measurements in this cycle');

ncwrite(filename,'/spectrometer1/last_sky_time',[calibratedSpectra.lastSkyTime]);
ncwriteatt(filename,'/spectrometer1/last_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/last_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/last_sky_time','description','stop time of the first sky measurements in this cycle');

ncwrite(filename,'/spectrometer1/time_min',[calibratedSpectra.timeMin]);
ncwriteatt(filename,'/spectrometer1/last_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/last_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/last_sky_time','description','minimum theoretical start time for this calibration cycle');

%%%%%%%%%%%%%%%%%
% Calibration variables
ncwrite(filename,'/spectrometer1/Tb',vertcat(calibratedSpectra.Tb)');
ncwrite(filename,'/spectrometer1/stdTb',vertcat(calibratedSpectra.stdTb)');
ncwrite(filename,'/spectrometer1/mean_std_Tb',[calibratedSpectra.meanStdTb]);
ncwrite(filename,'/spectrometer1/frequencies',calibratedSpectra(1).freq);

ncwrite(filename,'/spectrometer1/THot',[calibratedSpectra.THot]);
ncwrite(filename,'/spectrometer1/TSys',[calibratedSpectra.TSys]);
ncwrite(filename,'/spectrometer1/stdTSys',[calibratedSpectra.stdTSys]);
ncwrite(filename,'/spectrometer1/calibration_time',60*[calibratedSpectra.calibrationTime]);
ncwrite(filename,'/spectrometer1/mean_sky_elevation_angle',[calibratedSpectra.meanAngleAntenna]);

if isfield(calibratedSpectra,'if')
    ncwrite(filename,'/spectrometer1/intermediate_freq',calibratedSpectra(1).if);
else
    ncwrite(filename,'/spectrometer1/intermediate_freq',vertcat(-9999*ones(length(calibratedSpectra(1).freq),1))');
end

if isfield(calibratedSpectra,'stdTHot')
    ncwrite(filename,'/spectrometer1/stdTHot',[calibratedSpectra.stdTHot]);
else
    ncwrite(filename,'/spectrometer1/stdTHot',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'TempRoom')
    ncwrite(filename,'/spectrometer1/TRoom',[calibratedSpectra.TempRoom]);
    ncwrite(filename,'/spectrometer1/stdTRoom',[calibratedSpectra.stdTempRoom]);
else
    ncwrite(filename,'/spectrometer1/TRoom',-9999*ones(length(calibratedSpectra),1));
    ncwrite(filename,'/spectrometer1/stdTRoom',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'TempOut')
    ncwrite(filename,'/spectrometer1/TOut',[calibratedSpectra.TempOut]);
else
    ncwrite(filename,'/spectrometer1/TOut',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'TempWindow')
    ncwrite(filename,'/spectrometer1/TWindow',[calibratedSpectra.TempWindow]);
else
    ncwrite(filename,'/spectrometer1/TWindow',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'noiseLevel')
    ncwrite(filename,'/spectrometer1/noise_level',[calibratedSpectra.noiseLevel]);
else
    ncwrite(filename,'/spectrometer1/noise_level',-9999*ones(length(calibratedSpectra),1));
end

numInd=vertcat(calibratedSpectra.numberOfIndices);
ncwrite(filename,'/spectrometer1/number_of_hot_spectra',numInd(:,1));
ncwrite(filename,'/spectrometer1/number_of_cold_spectra',numInd(:,2));
ncwrite(filename,'/spectrometer1/number_of_sky_spectra',numInd(:,3));

ncwrite(filename,'/spectrometer1/mean_hot_counts',[calibratedSpectra.meanHotCounts]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the flags variables 
ncwrite(filename,'/flags/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/flags/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/flags/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/flags/time','description','mean time of the measurements for this cycle');

if isfield(calibratedSpectra,'errorVector')
    ncwrite(filename,'/flags/flags',1:length(calibratedSpectra(1).errorVector));
    ncwrite(filename,'/flags/calibration_flags',vertcat(calibratedSpectra.errorVector)');
else 
    ncwrite(filename,'/flags/flags',1);
    ncwrite(filename,'/flags/calibration_flags',-9999);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing Meteo Data
if isfield(logFile,'meteo') && ~isempty(fieldnames(logFile.meteo))
    ncwrite(filename,'/meteo/time',[logFile.meteo.dateNum]);
    ncwrite(filename,'/meteo/air_pressure',[logFile.meteo.air_pressure]);
    ncwrite(filename,'/meteo/air_temperature',[logFile.meteo.air_temperature]);
    ncwrite(filename,'/meteo/relative_humidity',[logFile.meteo.rel_humidity]);
    ncwrite(filename,'/meteo/precipitation',[logFile.meteo.precipitation]);
else
    ncwrite(filename,'/meteo/time',-9999);
    ncwrite(filename,'/meteo/air_pressure',-9999);
    ncwrite(filename,'/meteo/air_temperature',-9999);
    ncwrite(filename,'/meteo/relative_humidity',-9999);
    ncwrite(filename,'/meteo/precipitation',-9999);
end

% add time attributes directly
ncwriteatt(filename,'/meteo/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/meteo/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/meteo/time','description','time from the meteo stations');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originator attributes
% TODO --> CF convention ?
ncwriteatt(filename,'/','title','Brightness temperature measured by ground-based radiometer');
ncwriteatt(filename,'/','location',calibrationTool.dataLocation);
ncwriteatt(filename,'/','source',calibrationTool.dataSource);
ncwriteatt(filename,'/','name',calibrationTool.PI_NAME);
ncwriteatt(filename,'/','institution',calibrationTool.PI_AFFILIATION);
ncwriteatt(filename,'/','contact',calibrationTool.PI_ADDRESS);
ncwriteatt(filename,'/','mail',calibrationTool.PI_EMAIL);
ncwriteatt(filename,'/','instrument',calibrationTool.instrumentName);
ncwriteatt(filename,'/','number_of_spectrometer',calibrationTool.numberOfSpectrometer);
ncwriteatt(filename,'/','history','');
ncwriteatt(filename,'/','references','');
ncwriteatt(filename,'/','comment','');
ncwriteatt(filename,'/','raw_data',logFile.file);
if isfield(logFile,'SW_version')
    ncwriteatt(filename,'/','raw_data_software_version',num2str(logFile.SW_version(1)));
else
    ncwriteatt(filename,'/','raw_data_software_version','unknown');
end
ncwriteatt(filename,'/','calibration_version',calibratedSpectra(1).calibrationVersion);
ncwriteatt(filename,'/','raw_file_comment',logFile.comment);

ncwriteatt(filename,'/','raw_file_warning',warningLevel0);

if ~isempty(fieldnames(calibrationTool.labviewLog))
    if sum(isbetween([calibrationTool.labviewLog.dateTime], calibratedSpectra(1).theoreticalStartTime, calibratedSpectra(end).theoreticalStartTime + minutes(10))) > 0
        ncwriteatt(filename,'/','labview_logfile_warning','check labview log !');
    else
        ncwriteatt(filename,'/','labview_logfile_warning','clean');
    end
else
    ncwriteatt(filename,'/','labview_logfile_warning','no labview log found');
end

% Geolocation attributes
ncwriteatt(filename,'/','data_start_date',calibratedSpectra(1).firstSkyTime);
ncwriteatt(filename,'/','data_stop_date',calibratedSpectra(end).lastSkyTime);

ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(datetime('now','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','raw_filename','')
ncwriteatt(filename,'/','featureType','timeSeries');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attributes for groups and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the attributes of flags group

if isfield(calibratedSpectra,'errorVector')
    ncwriteatt(filename,'/flags','description','Each spectra is associated with a flag vector. The order and meaning of the flags are described in its attributes');
    ncwriteatt(filename,'/flags','number_of_flags',lenErrorVect);
    %ncwriteatt(filename,'/flags/calibration_flags','errorCode_1',calibratedSpectra(1).errorVectorDescription(1));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_2',calibratedSpectra(1).errorVectorDescription(2));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_3',calibratedSpectra(1).errorVectorDescription(3));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_4',calibratedSpectra(1).errorVectorDescription(4));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_5',calibratedSpectra(1).errorVectorDescription(5));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_6',calibratedSpectra(1).errorVectorDescription(6));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_7',calibratedSpectra(1).errorVectorDescription(7));
    for i = 1:lenErrorVect
        varName = ['errorCode_' num2str(i)];
        ncwriteatt(filename,'/flags/calibration_flags',varName,calibratedSpectra(1).errorVectorDescription(i));
    end
else
    ncwriteatt(filename,'/flags','description','No flags has been saved for this day...');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables attributes for the spectrometers and meteo group
% The following are required for each variable (CF convention):
ncwriteatt(filename,'/spectrometer1','spectrometer_type',calibrationTool.spectrometer);

attrName={'long_name','standard_name','units','description'};

attrVal.tod = {'TOD',...
    '',...
    'hour',...
    'Time of the day'};

attrVal.lat = {'station latitude',...
    'latitude',...
    'degree_north',...
    'latitude defined according to WGS84'};

attrVal.lon = {'station longitude',...
    'longitude',...
    'degree_east',...
    'longitude defined according to WGS84'};

attrVal.alt = {'station altitude',...
    'altitude',...
    'm',...
    'above see level'};

attrVal.azimuth = {'azimuth angle',...
    'sensor_azimuth_angle',...
    'degree',...
    'angle measured clockwise positive, 0 deg is northwise'};

%attrVal.effCalTime = {'effective calibration time',...
%    'calibration_time',...
%    'second',...
 %   'still to improve'};

attrVal.Tb = {'Tb',...
    'brightness_temperature',...
    'K',...
    'calibrated brightness temperature for this cycle'};

attrVal.stdTb = {'stdTb',...
    'spectra of stdTb',...
    'K',...
    'standard deviation of brightness temperature for this cycle per channel'};

attrVal.freq = {'f',...
    'frequency vector',...
    'Hz',...
    'frequency vector for the spectrometer'};

attrVal.meanStdTb = {'mean_std_Tb',...
    'mean stdTb',...
    'K',...
    'mean standard deviation of brightness temperature for this cycle (without bad channel)'};

attrVal.if = {'if',...
    'intermediate frequency vector',...
    'Hz',...
    'intermediate frequency vector for the spectrometer'};

attrVal.THot = {'THot',...
    '',...
    'K',...
    'Mean temperature of the hot load'};

attrVal.stdTHot = {'stdTHot',...
    '',...
    'K',...
    'standard deviation of the hot load temperature'};

attrVal.TSys = {'TSys',...
    '',...
    'K',...
    'mean system temperature'};

attrVal.stdTSys = {'stdTSys',...
    'standard deviation of the system temperature',...
    'K',...
    ''};

attrVal.calibrationTime = {'calibrationTime',...
    '',...
    'second',...
    'Time used for calibrating the spectra'};

attrVal.noiseLevel = {'noise_level',...
    'std(diff(Tb))',...
    'Tb',...
    'describes how noisy is the spectra'};

attrVal.meanAngleAntenna = {'meanAngleAntenna',...
    'elevation_angle',...
    'degree',...
    'mean elevation angle of the antenna during this cycle'};

attrVal.TRoom = {'TRoom',...
    '',...
    'K',...
    'mean room temperature'};

attrVal.stdTRoom = {'stdTRoom',...
    '',...
    'K',...
    'standard deviation of room temperature'};

attrVal.TOut = {'TOut',...
    '',...
    'K',...
    'mean outside temperature'};

attrVal.TWindow = {'TWindow',...
    '',...
    'K',...
    'mean window temperature'};

% for Meteo data:
attrVal.airP = {'air_pressure',...
    '',...
    'hPa',...
    'air pressure at the station'};

attrVal.airT = {'air_temperature',...
    '',...
    'K',...
    'air temperature at the station'};

attrVal.relH = {'relative_humidity',...
    '',...
    '1',...
    'relative humidity of the air at the station'};

attrVal.precipitation = {'precipitation',...
    '',...
    '...',...
    'precipitation'};

% Ugly and open to suggestion
for i=1:length(attrName)
    ncwriteatt(filename,'/spectrometer1/time_of_day',attrName{i},attrVal.tod{i});
    ncwriteatt(filename,'/spectrometer1/lat',attrName{i},attrVal.lat{i});
    ncwriteatt(filename,'/spectrometer1/lon',attrName{i},attrVal.lon{i});
    ncwriteatt(filename,'/spectrometer1/alt',attrName{i},attrVal.alt{i});
    ncwriteatt(filename,'/spectrometer1/azimuth_angle',attrName{i},attrVal.azimuth{i});
    %ncwriteatt(filename,'/spectrometer1/effectiveCalibrationTime',attrName{i},attrVal.effCalTime{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/stdTb',attrName{i},attrVal.stdTb{i});
    ncwriteatt(filename,'/spectrometer1/mean_std_Tb',attrName{i},attrVal.meanStdTb{i});
    ncwriteatt(filename,'/spectrometer1/frequencies',attrName{i},attrVal.freq{i});
    ncwriteatt(filename,'/spectrometer1/intermediate_freq',attrName{i},attrVal.if{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/TSys',attrName{i},attrVal.TSys{i});
    ncwriteatt(filename,'/spectrometer1/stdTSys',attrName{i},attrVal.stdTSys{i});
    ncwriteatt(filename,'/spectrometer1/calibration_time',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/mean_sky_elevation_angle',attrName{i},attrVal.meanAngleAntenna{i});
    ncwriteatt(filename,'/spectrometer1/TRoom',attrName{i},attrVal.TRoom{i});
    ncwriteatt(filename,'/spectrometer1/stdTRoom',attrName{i},attrVal.stdTRoom{i});
    ncwriteatt(filename,'/spectrometer1/TWindow',attrName{i},attrVal.TWindow{i});
    ncwriteatt(filename,'/spectrometer1/TOut',attrName{i},attrVal.TOut{i});
    ncwriteatt(filename,'/spectrometer1/noise_level',attrName{i},attrVal.noiseLevel{i});
    
    % Meteo attr
    ncwriteatt(filename,'/meteo/air_pressure',attrName{i},attrVal.airP{i});
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.airT{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relH{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding debug group (optionnal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are 2 types of debug information for which we create new groups.
if calibrationTool.calType=="debug"
    nccreate(filename,'/UpDown/timeUpDown','Dimensions',{'timeUpDown',Inf},'Datatype','double','Format','netcdf4');
    nccreate(filename,'/UpDown/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','int64','FillValue',-9999)
    
    nccreate(filename,'/UpDown/TbUp','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'timeUpDown',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/UpDown/TbDown','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'timeUpDown',Inf},'Datatype','double','FillValue',-9999);

    ncwrite(filename,'/UpDown/TbUp',vertcat(calibratedSpectra.TbUp)');
    ncwrite(filename,'/UpDown/TbDown',vertcat(calibratedSpectra.TbDown)');
    
    % For all cycle:
    % initialize matrices
    nClean=0;
    for t = 1:length(calibratedSpectra)
        nClean=nClean+calibratedSpectra(t).numberOfCleanAntennaAngle;
    end
    
    channel_idx=int64(ones(nClean,calibrationTool.numberOfChannels)*NaN);
    TbAll=single(ones(nClean,calibrationTool.numberOfChannels)*NaN);
    cycleNumber=int64(ones(nClean,1)*NaN);
    cycleId=int64(1:nClean);
    
    counter=1;
    for t = 1:length(calibratedSpectra)
        for i = 1:calibratedSpectra(t).numberOfCleanAntennaAngle
            TbAll(counter,:)=calibratedSpectra(t).TbAll(i,:);
            channel_idx(counter,:)=1:calibrationTool.numberOfChannels;
            cycleNumber(counter)=t;
            counter=counter+1;
        end
    end
    
    % Group for debugging variables:
    nccreate(filename,'/debug/cycle_id','Dimensions',{'cycle_id',Inf},'Datatype','int64')
    nccreate(filename,'/debug/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','int64','FillValue',-9999)
    
    nccreate(filename,'/debug/calibration_cycle_number','Dimensions',{'cycle_id',Inf},'Datatype','int64','FillValue',-9999)
    nccreate(filename,'/debug/Tb_all','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'cycle_id',Inf},'Datatype','single','FillValue',-9999)
    
    % Writing the debug variables
    ncwrite(filename,'/debug/cycle_id',cycleId);
    
    ncwrite(filename,'/debug/channel_idx',1:calibrationTool.numberOfChannels);
    ncwrite(filename,'/debug/calibration_cycle_number',cycleNumber');
    ncwrite(filename,'/debug/Tb_all',TbAll');   
end

disp(['File saved as: ' filename])
end
