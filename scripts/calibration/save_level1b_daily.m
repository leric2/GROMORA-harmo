function calibrationTool = save_level1b_daily(calibrationTool,level1b)
%==========================================================================
% NAME          | save_level1a_daily.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Function saving the integrated spectra for a complete day
%               | in the form of netCDF4 file. 
%               | The file contains one entry
%               | per integration cylce (typically 1 hour) on a time 
%               | coordinate as well as all necessery metadata to identify
%               | and re-use the data.
%               |
%               |
% ARGUMENTS     | INPUTS: - calibrationTool
%               |         - logFile
%               |         - calibratedSpectra
%               |         - warningLevel0
%               |         
%               | OUTPUTS: - calibrationTool with added field
%               |          
%               |
% SAVE          | Level1b netCDF4 containing the following groups:
%               |   - '/' for global attributes
%               |   - '/spectrometer1/' for integration
%               |   - '/flags/' for
%               |   - '/meteo/' 
%               |
%               |
%==========================================================================

% Filename and location for DAILY netCDF file
locationLevel1b=calibrationTool.level1Folder;
integratedSpectra=level1b.integration;

if calibrationTool.integrationTime == 60
    filename=[locationLevel1b calibrationTool.instrumentName '_level1b_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];
else
    filename=[locationLevel1b calibrationTool.instrumentName '_level1b_' num2str(calibrationTool.integrationTime/60) 'h_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];
end
calibrationTool.filenameLevel1b=filename;

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
nccreate(filename,'/spectrometer1/year','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/month','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/day','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/time_of_day','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/first_sky_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/last_sky_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/time_min','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);

%%%%%%%%%%%%%%%%%
% Integration variables   
%nccreate(filename,'/spectrometer1/effectiveCalibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/Tb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/Tb_corr','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/stdTb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/good_channels','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);

nccreate(filename,'/spectrometer1/frequencies','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/intermediate_freq','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/integration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_sky_elevation_angle','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_std_Tb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/TRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TWindow','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_calibrated_spectra','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_of_hot_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_cold_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_sky_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)

% Tropospheric correction data:
nccreate(filename,'/spectrometer1/tropospheric_transmittance','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/tropospheric_opacity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
nccreate(filename,'/meteo/time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/air_pressure','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/air_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/relative_humidity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/precipitation','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags dataset

% Coordinates time and flags
nccreate(filename,'/flags/time','Dimensions',{'time',Inf},'Datatype','double')
if isfield(integratedSpectra,'errorVector')
    lenErrorVect=length(integratedSpectra(1).errorVector);
else
    lenErrorVect=1;
end

nccreate(filename,'/flags/flags','Dimensions',{'flags',lenErrorVect},'Datatype','int64')
nccreate(filename,'/flags/calibration_flags','Dimensions',{'flags',lenErrorVect,'time',Inf},'Datatype','int64','FillValue',-9999)
%nccreate(filename,'/error/errorLabel','Dimensions',{'errorLabel',20,'time',Inf},'Datatype','char','FillValue','')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the spectrometer1 variable
% Coordinate variables, directly adding the attributes
ncwrite(filename,'/spectrometer1/time',[integratedSpectra.meanDatetime]);
ncwriteatt(filename,'/spectrometer1/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/time','description','mean time of THE BEGINNING of all antenna measurements for this cycle');

ncwrite(filename,'/spectrometer1/channel_idx',1:calibrationTool.numberOfChannels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some variables to help geolocate the file:
% these are scalar variables as they do not vary in time
ncwrite(filename,'/spectrometer1/lat',ones(length(integratedSpectra),1)*calibrationTool.lat);
ncwrite(filename,'/spectrometer1/lon',ones(length(integratedSpectra),1)*calibrationTool.lon);
ncwrite(filename,'/spectrometer1/alt',ones(length(integratedSpectra),1)*calibrationTool.altitude);
ncwrite(filename,'/spectrometer1/azimuth_angle',ones(length(integratedSpectra),1)*calibrationTool.azimuthAngle);

% some variable for better identifying the time period of the measurements
ncwrite(filename,'/spectrometer1/year',int64([integratedSpectra.year]));
ncwrite(filename,'/spectrometer1/month',int64([integratedSpectra.month]));
ncwrite(filename,'/spectrometer1/day',int64([integratedSpectra.day]));
ncwrite(filename,'/spectrometer1/time_of_day',[integratedSpectra.TOD]);

ncwrite(filename,'/spectrometer1/first_sky_time',[integratedSpectra.firstSkyTime]);  
ncwriteatt(filename,'/spectrometer1/first_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/first_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/first_sky_time','description','start time of the first sky measurements in this cycle');

ncwrite(filename,'/spectrometer1/last_sky_time',[integratedSpectra.lastSkyTime]);
ncwriteatt(filename,'/spectrometer1/last_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/last_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/last_sky_time','description','stop time of the first sky measurements in this cycle');

ncwrite(filename,'/spectrometer1/time_min',[integratedSpectra.timeMin]);
ncwriteatt(filename,'/spectrometer1/last_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/last_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/last_sky_time','description','minimum theoretical start time for this calibration cycle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main set of variables
ncwrite(filename,'/spectrometer1/Tb',vertcat(integratedSpectra.Tb)');
ncwrite(filename,'/spectrometer1/Tb_corr',vertcat(integratedSpectra.TbTroposphericWindowCorr)');
ncwrite(filename,'/spectrometer1/stdTb',vertcat(integratedSpectra.stdTb)');
ncwrite(filename,'/spectrometer1/good_channels',~vertcat(integratedSpectra.potentialBadChannels)');
ncwrite(filename,'/spectrometer1/frequencies',integratedSpectra(1).freq);
ncwrite(filename,'/spectrometer1/intermediate_freq',integratedSpectra(1).if);
ncwrite(filename,'/spectrometer1/THot',[integratedSpectra.THot]);
%ncwrite(filename,'/spectrometer1/stdTHot',[integratedSpectra.stdTHot]);
ncwrite(filename,'/spectrometer1/TSys',[integratedSpectra.TSys]);
%ncwrite(filename,'/spectrometer1/stdTSys',[integratedSpectra.stdTSys]);
ncwrite(filename,'/spectrometer1/calibration_time',[integratedSpectra.calibrationTime]);
ncwrite(filename,'/spectrometer1/integration_time',[integratedSpectra.integrationTime]);
ncwrite(filename,'/spectrometer1/mean_sky_elevation_angle',[integratedSpectra.meanAngleAntenna]);
ncwrite(filename,'/spectrometer1/mean_std_Tb',[integratedSpectra.meanStdTbFromCal]);

ncwrite(filename,'/spectrometer1/number_calibrated_spectra',[integratedSpectra.numberOfAveragedSpectra]);

% Tropospheric correction data:
ncwrite(filename,'/spectrometer1/tropospheric_transmittance',[integratedSpectra.troposphericTransmittance]);
ncwriteatt(filename,'/spectrometer1/tropospheric_transmittance','method',integratedSpectra(1).troposphericCorrType);

ncwrite(filename,'/spectrometer1/tropospheric_opacity',[integratedSpectra.troposphericOpacity]);
ncwriteatt(filename,'/spectrometer1/tropospheric_opacity','method',integratedSpectra(1).troposphericCorrType);

% Data that are not present for every instrument
if isfield(integratedSpectra,'TempRoom')
    ncwrite(filename,'/spectrometer1/TRoom',[integratedSpectra.TempRoom]);
    ncwrite(filename,'/spectrometer1/stdTRoom',[integratedSpectra.stdTempRoom]);
else
    ncwrite(filename,'/spectrometer1/TRoom',-9999*ones(length(integratedSpectra),1));
    ncwrite(filename,'/spectrometer1/stdTRoom',-9999*ones(length(integratedSpectra),1));
end

if isfield(integratedSpectra,'TOut')
    ncwrite(filename,'/spectrometer1/TOut',[integratedSpectra.TOut]);
else
    ncwrite(filename,'/spectrometer1/TOut',-9999*ones(length(integratedSpectra),1));
end

if isfield(integratedSpectra,'TWindow')
    ncwrite(filename,'/spectrometer1/TWindow',[integratedSpectra.TWindow]);
else
    ncwrite(filename,'/spectrometer1/TWindow',-9999*ones(length(integratedSpectra),1));
end

ncwrite(filename,'/spectrometer1/number_of_hot_spectra',[integratedSpectra.numHotSpectra]);
ncwrite(filename,'/spectrometer1/number_of_cold_spectra',[integratedSpectra.numColdSpectra]);
ncwrite(filename,'/spectrometer1/number_of_sky_spectra',[integratedSpectra.numSkySpectra]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
ncwrite(filename,'/meteo/time',[integratedSpectra.meanDatetime]);
ncwriteatt(filename,'/meteo/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/meteo/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/meteo/time','description','mean time of THE BEGINNING of all antenna measurements for this cycle');

ncwrite(filename,'/meteo/air_temperature',[integratedSpectra.meanAirTemperature]);
ncwrite(filename,'/meteo/air_pressure',[integratedSpectra.meanAirPressure]);
ncwrite(filename,'/meteo/relative_humidity',[integratedSpectra.meanRelativeHumidity]);
ncwrite(filename,'/meteo/precipitation',[integratedSpectra.rainAccumulation]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the flags variables

ncwrite(filename,'/flags/time',[integratedSpectra.meanDatetime]);
ncwriteatt(filename,'/flags/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/flags/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/flags/time','description','mean time of the measurements for this cycle');

if isfield(integratedSpectra,'errorVector')
    ncwrite(filename,'/flags/flags',1:length(integratedSpectra(1).errorVector));
    ncwrite(filename,'/flags/calibration_flags',vertcat(integratedSpectra.errorVector)');
else 
    ncwrite(filename,'/flags/flags',1);
    ncwrite(filename,'/flags/calibration_flags',-9999);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originator attributes
ncwriteatt(filename,'/','title','Brightness temperature measured by ground-based radiometer');
ncwriteatt(filename,'/','location',calibrationTool.dataLocation);
ncwriteatt(filename,'/','source',calibrationTool.dataSource);
ncwriteatt(filename,'/','name',calibrationTool.PI_NAME);
ncwriteatt(filename,'/','institution',calibrationTool.PI_AFFILIATION);
ncwriteatt(filename,'/','contact',calibrationTool.PI_ADDRESS);
ncwriteatt(filename,'/','mail',calibrationTool.PI_EMAIL);
ncwriteatt(filename,'/','history','');
ncwriteatt(filename,'/','references','');
ncwriteatt(filename,'/','comment','');
%ncwriteatt(filename,'/','DATA_VARIABLES','');
ncwriteatt(filename,'/','instrument',calibrationTool.instrumentName);
ncwriteatt(filename,'/','number_of_spectrometer',calibrationTool.numberOfSpectrometer);
ncwriteatt(filename,'/','raw_data_filename',calibrationTool.logFile.rawFilename);
ncwriteatt(filename,'/','raw_data_software_version',calibrationTool.logFile.raw_data_software_version);
ncwriteatt(filename,'/','calibration_version',calibrationTool.logFile.calibration_version);
ncwriteatt(filename,'/','raw_file_comment',calibrationTool.logFile.raw_file_comment);
ncwriteatt(filename,'/','comment',calibrationTool.logFile.comment);
ncwriteatt(filename,'/','filename_level1a',calibrationTool.logFile.filenameLevel1a);
ncwriteatt(filename,'/','creation_date_level1a',calibrationTool.logFile.creation_date_level1a);
ncwriteatt(filename,'/','raw_file_warning',calibrationTool.logFile.raw_file_warning);

% Global file attributes
ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(now,'yyyymmddTHHMMSSZ'));   %TODO
ncwriteatt(filename,'/','featureType','timeSeries');

if ~isempty(fieldnames(calibrationTool.labviewLog))
    if sum(isbetween([calibrationTool.labviewLog.dateTime],...
            datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day),...
            datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day+1))) > 0
        ncwriteatt(filename,'/','labview_logfile_warning','check labview log !');
    else
        ncwriteatt(filename,'/','labview_logfile_warning','clean');
    end
else
    ncwriteatt(filename,'/','labview_logfile_warning','no labview log found');
end

% Geolocation attributes
%ncwriteatt(filename,'/','data_start_date',datestr(calibratedSpectra(1).dateStart,'yyyymmddTHHMMSSZ'));
%ncwriteatt(filename,'/','data_stop_date',datestr(calibratedSpectra(end).dateStop,'yyyymmddTHHMMSSZ'));

%ncwriteatt(filename,'/','data_start_date',integratedSpectra(1).dateStart);
%ncwriteatt(filename,'/','data_stop_date',integratedSpectra(end).dateStop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attributes for groups and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the attributes of flags group
if isfield(integratedSpectra,'errorVector')
    ncwriteatt(filename,'/flags','description','Each spectra is associated with a flag vector. The order and meaning of the flags are described in its attributes');
    for i = 1:lenErrorVect
        varName = ['errorCode_' num2str(i)];
        ncwriteatt(filename,'/flags/calibration_flags',varName,integratedSpectra(1).errorVectorDescription(i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables attributes for the spectrometers group
% The following are required for each variable
% Attribute name (CF convention)
attrName={'long_name','standard_name','units','description'};

%Attributes for the spectrometer1 variables (CF convention)
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
    'integrated brightness temperature for this cycle'};

attrVal.Tb_corr = {'Tb_corr',...
    'corrected brightness_temperature',...
    'K',...
    'integrated brightness temperature for this cycle, corrected for window and troposphere'};

attrVal.stdTb = {'stdTb',...
    'spectra of stdTb',...
    'K',...
    'standard deviation of brightness temperature for this cycle per channel'};

attrVal.good_channels = {'good_channels',...
    'good_channels',...
    '1',...
    'Flags for the quality of the channels performed during the integration'};

attrVal.freq = {'f',...
    'frequency vector',...
    'Hz',...
    'frequency vector for the spectrometer'};

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

attrVal.meanStdTb = {'mean_std_Tb',...
    'mean stdTb',...
    'K',...
    'mean standard deviation of brightness temperature for this integration cycle (without bad channel)'};

attrVal.calibrationTime = {'calibrationTime',...
    '',...
    'second',...
    'Time used for calibrating the spectra'};

attrVal.integrationTime = {'integrationTime',...
    '',...
    'second',...
    'Time used for integrating the spectra'};

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

% Tropospheric correction data:
attrVal.tropospheric_transmittance = {'tropospheric transmittance',...
    '',...
    '1',...
    ''};

attrVal.tropospheric_opacity = {'tropospheric opacity',...
    '',...
    '1',...
    ''};

% Meteo Data
% Tropospheric correction data:
attrVal.air_temperature = {'air temperature',...
    'air_temperature',...
    'K',...
    ''};

attrVal.air_pressure = {'air pressure',...
    'air_pressure',...
    'hPa',...
    'Pressure at the station'};

attrVal.relative_humidity = {'relative humidity',...
    'relative_humidity',...
    '1',...
    ''};

attrVal.precipitation = {'rain accumulation',...
    '',...
    'mm',...
    ''};

for i=1:length(attrName)
    ncwriteatt(filename,'/spectrometer1/time_of_day',attrName{i},attrVal.tod{i});
    ncwriteatt(filename,'/spectrometer1/lat',attrName{i},attrVal.lat{i});
    ncwriteatt(filename,'/spectrometer1/lon',attrName{i},attrVal.lon{i});
    ncwriteatt(filename,'/spectrometer1/alt',attrName{i},attrVal.alt{i});
    ncwriteatt(filename,'/spectrometer1/azimuth_angle',attrName{i},attrVal.azimuth{i});
    %ncwriteatt(filename,'/spectrometer1/effectiveCalibrationTime',attrName{i},attrVal.effCalTime{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/Tb_corr',attrName{i},attrVal.Tb_corr{i});
    ncwriteatt(filename,'/spectrometer1/stdTb',attrName{i},attrVal.stdTb{i});
    ncwriteatt(filename,'/spectrometer1/good_channels',attrName{i},attrVal.good_channels{i});
    ncwriteatt(filename,'/spectrometer1/frequencies',attrName{i},attrVal.freq{i});
    ncwriteatt(filename,'/spectrometer1/intermediate_freq',attrName{i},attrVal.if{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    %ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/TSys',attrName{i},attrVal.TSys{i});
    ncwriteatt(filename,'/spectrometer1/mean_std_Tb',attrName{i},attrVal.meanStdTb{i});
    ncwriteatt(filename,'/spectrometer1/calibration_time',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/integration_time',attrName{i},attrVal.integrationTime{i});
    ncwriteatt(filename,'/spectrometer1/mean_sky_elevation_angle',attrName{i},attrVal.meanAngleAntenna{i});
    ncwriteatt(filename,'/spectrometer1/TRoom',attrName{i},attrVal.TRoom{i});
    ncwriteatt(filename,'/spectrometer1/stdTRoom',attrName{i},attrVal.stdTRoom{i});
    ncwriteatt(filename,'/spectrometer1/TWindow',attrName{i},attrVal.TWindow{i});
    ncwriteatt(filename,'/spectrometer1/TOut',attrName{i},attrVal.TOut{i});
    
    % Meteo
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.air_temperature{i});
    ncwriteatt(filename,'/meteo/air_pressure',attrName{i},attrVal.air_pressure{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relative_humidity{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});
    
    % Corrections
    ncwriteatt(filename,'/spectrometer1/tropospheric_transmittance',attrName{i},attrVal.tropospheric_transmittance{i});
    ncwriteatt(filename,'/spectrometer1/tropospheric_opacity',attrName{i},attrVal.tropospheric_opacity{i});
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Adding debug (optionnal)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if calibrationTool.calType=="debug"
%     % initialize matrices
%     nClean=0;
%     for t = 1:length(integratedSpectra)
%         nClean=nClean+integratedSpectra(t).numberOfCleanAntennaAngle;
%     end
%     
%     channel_idx=int64(ones(nClean,calibrationTool.numberOfChannels)*NaN);
%     TbAll=single(ones(nClean,calibrationTool.numberOfChannels)*NaN);
%     cycleNumber=int64(ones(nClean,1)*NaN);
%     cycleId=int64(1:nClean);
%     
%     counter=1;
%     for t = 1:length(integratedSpectra)
%         for i =1:integratedSpectra(t).numberOfCleanAntennaAngle
%             TbAll(counter,:)=integratedSpectra(t).TbAll(i,:);
%             channel_idx(counter,:)=1:calibrationTool.numberOfChannels;
%             cycleNumber(counter)=t;
%             counter=counter+1;
%         end
%     end
%     
%     % Group for debugging variables:
%     nccreate(filename,'/debug/cycle_id','Dimensions',{'cycle_id',Inf},'Datatype','int64')
%     nccreate(filename,'/debug/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','int64','FillValue',-9999)
%     
%     nccreate(filename,'/debug/calibration_cycle_number','Dimensions',{'cycle_id',Inf},'Datatype','int64','FillValue',-9999)
%     nccreate(filename,'/debug/Tb_all','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'cycle_id',Inf},'Datatype','single','FillValue',-9999)
%     
%     % Writing the debug variables
%     ncwrite(filename,'/debug/cycle_id',cycleId);
%     
%     ncwrite(filename,'/debug/channel_idx',1:calibrationTool.numberOfChannels);
%     ncwrite(filename,'/debug/calibration_cycle_number',cycleNumber');
%     ncwrite(filename,'/debug/Tb_all',TbAll');
% end

disp(['File saved as: ' filename])
end
