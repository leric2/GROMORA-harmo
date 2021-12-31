function calibrationTool = save_level1b_daily(calibrationTool,level1)
%==========================================================================
% NAME      | save_level1b_daily.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 01.2020
%           |
% ABSTRACT  | Function saving the integrated spectra for a complete day
%           | in the form of netCDF4 file. The file contains one entry 
%           | per integration cylce (typically 1 hour) on a time 
%           | coordinate as well as all necessery metadata to identify
%           | and re-use the data.
%           |
%           |
% ARGUMENTS | INPUTS: 1. calibrationTool:
%           |           - level1Folder
%           |           - instrumentName
%           |           - spectrometer
%           |           - dateStr
%           |           - integrationTime
%           |           - numberOfChannels
%           |           - meanDatetimeUnit
%           |           - calendar
%           |           - lat, lon, altitude, azimuthAngle
%           |           - dataLocation, dataSource
%           |           - PI_NAME, PI_AFFILIATION, PI_ADDRESS, PI_EMAIL
%           |           - numberOfSpectrometer
%           |           - logFile: substructure containing metadata
%           |           - labviewLog, Day, Month, Year  (optional)
%           |         2. level1: structure containing both level 1a and
%           |            level 1b for this day.
%           | 
%           | OUTPUTS: - calibrationTool with added field filenameLevel1b
%           |      
% SAVE      | Level1b netCDF4 containing the following groups:
%           |   - '/' for global attributes
%           |   - '/spectrometer1/' for spectrometer dataset
%           |   - '/flags/'
%           |   - '/meteo/'
%           |  
%==========================================================================
% Filename and location for DAILY netCDF file
locationLevel1b=calibrationTool.level1Folder;
integratedSpectra=level1.integratedSpectra;

if calibrationTool.integrationTime == 60
    filename=[locationLevel1b calibrationTool.instrumentName '_level1b_' calibrationTool.spectrometer '_' calibrationTool.dateStr calibrationTool.extraName '.nc'];
else
    filename=[locationLevel1b calibrationTool.instrumentName '_level1b_' num2str(calibrationTool.integrationTime/60) 'h_' calibrationTool.spectrometer '_' calibrationTool.dateStr calibrationTool.extraName '.nc'];
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
nccreate(filename,'/spectrometer1/MJD2K','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
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
nccreate(filename,'/spectrometer1/Tb_win_corr','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/Tb_corr','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/stdTb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/good_channels','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);

if calibrationTool.savePlanckIntensity
    nccreate(filename,'/spectrometer1/intensity_planck','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/spectrometer1/intensity_planck_win_corr','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
end

nccreate(filename,'/spectrometer1/frequencies','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/intermediate_freq','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
%nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/noise_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
%nccreate(filename,'/spectrometer1/stdTNoise','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/integration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_sky_elevation_angle','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/meanTb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_std_Tb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/noise_level','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/TRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TWindow','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/VGunn','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_of_calibrated_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_of_hot_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_cold_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_sky_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)

% Tropospheric correction data:
nccreate(filename,'/spectrometer1/tropospheric_transmittance','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/tropospheric_opacity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/tropospheric_opacity_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

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
ncwriteatt(filename,'/spectrometer1/time','timezone',calibrationTool.timeZone);
ncwriteatt(filename,'/spectrometer1/time','description','mean time recorded at the beginning of all sky measurements during this integration cycle');
if isfield(integratedSpectra, 'MJD2K')
    ncwrite(filename,'/spectrometer1/MJD2K',[integratedSpectra.MJD2K]);
    ncwriteatt(filename,'/spectrometer1/MJD2K','units','MJD2K');
    ncwriteatt(filename,'/spectrometer1/MJD2K','calendar','Julian');
    ncwriteatt(filename,'/spectrometer1/MJD2K','description','mean time recorded at the beginning of all sky measurements during this integration cycle');
end
ncwrite(filename,'/spectrometer1/channel_idx',1:calibrationTool.numberOfChannels);
ncwriteatt(filename,'/spectrometer1/channel_idx','description','index of the spectrometer channels, from 1 to N (number of channels)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some variables to help geolocate the file:
% these are scalar variables as they do not vary in time
ncwrite(filename,'/spectrometer1/lat',ones(length(integratedSpectra),1)*calibrationTool.lat);
ncwrite(filename,'/spectrometer1/lon',ones(length(integratedSpectra),1)*calibrationTool.lon);
ncwrite(filename,'/spectrometer1/alt',ones(length(integratedSpectra),1)*calibrationTool.altitude);
ncwrite(filename,'/spectrometer1/azimuth_angle',ones(length(integratedSpectra),1)*calibrationTool.azimuthAngle);

% some variable for better identifying the time period of the measurements
ncwrite(filename,'/spectrometer1/year',int64([integratedSpectra.year]));
ncwriteatt(filename,'/spectrometer1/year','description','year of the measurement as integer');

ncwrite(filename,'/spectrometer1/month',int64([integratedSpectra.month]));
ncwriteatt(filename,'/spectrometer1/month','description','month of the measurement as integer');

ncwrite(filename,'/spectrometer1/day',int64([integratedSpectra.day]));
ncwriteatt(filename,'/spectrometer1/day','description','day of the month as integer');

ncwrite(filename,'/spectrometer1/time_of_day',[integratedSpectra.time_of_day]);

ncwrite(filename,'/spectrometer1/first_sky_time',[integratedSpectra.first_sky_time]);  
ncwriteatt(filename,'/spectrometer1/first_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/first_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/first_sky_time','description','time of the first sky measurements in this integration cycle');

ncwrite(filename,'/spectrometer1/last_sky_time',[integratedSpectra.last_sky_time]);
ncwriteatt(filename,'/spectrometer1/last_sky_time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/last_sky_time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/last_sky_time','description','time of the last sky measurements in this integration cycle');

ncwrite(filename,'/spectrometer1/time_min',[integratedSpectra.time_min]);
ncwriteatt(filename,'/spectrometer1/time_min','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time_min','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/time_min','description','minimum theoretical start time for this integration cycle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main set of variables
ncwrite(filename,'/spectrometer1/Tb',vertcat(integratedSpectra.Tb)');
ncwrite(filename,'/spectrometer1/Tb_win_corr',vertcat(integratedSpectra.TbWinCorr)');
ncwrite(filename,'/spectrometer1/Tb_corr',vertcat(integratedSpectra.TbTroposphericWindowCorr)');
ncwrite(filename,'/spectrometer1/stdTb',vertcat(integratedSpectra.stdTb)');
ncwrite(filename,'/spectrometer1/good_channels',~vertcat(integratedSpectra.potentialBadChannels)');
ncwrite(filename,'/spectrometer1/frequencies',vertcat(integratedSpectra.frequencies)');
ncwrite(filename,'/spectrometer1/intermediate_freq',vertcat(integratedSpectra.intermediate_freq)');
ncwrite(filename,'/spectrometer1/THot',[integratedSpectra.THot]);
%ncwrite(filename,'/spectrometer1/stdTHot',[integratedSpectra.stdTHot]);
ncwrite(filename,'/spectrometer1/noise_temperature',[integratedSpectra.TNoise]);
%ncwrite(filename,'/spectrometer1/stdTNoise',[integratedSpectra.stdTNoise]);
ncwrite(filename,'/spectrometer1/calibration_time',[integratedSpectra.calibration_time]);
ncwrite(filename,'/spectrometer1/integration_time',[integratedSpectra.integration_time]);
ncwrite(filename,'/spectrometer1/mean_sky_elevation_angle',[integratedSpectra.mean_sky_elevation_angle]);
ncwrite(filename,'/spectrometer1/mean_std_Tb',[integratedSpectra.meanStdTbFromCal]);
ncwrite(filename,'/spectrometer1/meanTb',[integratedSpectra.meanTb]);

ncwrite(filename,'/spectrometer1/number_of_calibrated_spectra',[integratedSpectra.numberOfAveragedSpectra]);

if calibrationTool.savePlanckIntensity
    ncwrite(filename,'/spectrometer1/intensity_planck',vertcat(integratedSpectra.intensity_planck)');
    ncwrite(filename,'/spectrometer1/intensity_planck_win_corr',vertcat(integratedSpectra.intensityPlanckWinCorr)');
end

% Tropospheric correction data:
ncwrite(filename,'/spectrometer1/tropospheric_transmittance',[integratedSpectra.troposphericTransmittance]);
ncwriteatt(filename,'/spectrometer1/tropospheric_transmittance','method',integratedSpectra(1).troposphericCorrType);

ncwrite(filename,'/spectrometer1/tropospheric_opacity',[integratedSpectra.troposphericOpacity]);
ncwriteatt(filename,'/spectrometer1/tropospheric_opacity','method',integratedSpectra(1).troposphericCorrType);

ncwrite(filename,'/spectrometer1/tropospheric_opacity_tc',[integratedSpectra.tropospheric_opacity_tc]);

if isfield(integratedSpectra,'noiseLevel')
    ncwrite(filename,'/spectrometer1/noise_level',[integratedSpectra.noiseLevel]);
else
    ncwrite(filename,'/spectrometer1/noise_level',-9999*ones(length(integratedSpectra),1));
end

% Data that are not present for every instrument
if isfield(integratedSpectra,'TRoom')
    ncwrite(filename,'/spectrometer1/TRoom',[integratedSpectra.TRoom]);
    ncwrite(filename,'/spectrometer1/stdTRoom',-9999*ones(length(integratedSpectra),1));
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

ncwrite(filename,'/spectrometer1/VGunn',[integratedSpectra.VGunn]);

ncwrite(filename,'/spectrometer1/number_of_hot_spectra',[integratedSpectra.number_of_hot_spectra]);
ncwrite(filename,'/spectrometer1/number_of_cold_spectra',[integratedSpectra.number_of_cold_spectra]);
ncwrite(filename,'/spectrometer1/number_of_sky_spectra',[integratedSpectra.number_of_sky_spectra]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
ncwrite(filename,'/meteo/time',[integratedSpectra.meanDatetime]);
ncwriteatt(filename,'/meteo/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/meteo/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/meteo/time','description','mean time recorded at the beginning of all sky measurements during this calibration cycle');
ncwriteatt(filename,'/meteo/time','timezone',calibrationTool.timeZone);

ncwrite(filename,'/meteo/air_temperature',[integratedSpectra.mean_air_temperature]);
ncwrite(filename,'/meteo/air_pressure',[integratedSpectra.mean_air_pressure]);
ncwrite(filename,'/meteo/relative_humidity',[integratedSpectra.mean_relative_humidity]);
ncwrite(filename,'/meteo/precipitation',[integratedSpectra.rain_accumulation]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the flags variables

ncwrite(filename,'/flags/time',[integratedSpectra.meanDatetime]);
ncwriteatt(filename,'/flags/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/flags/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/flags/time','description','mean time of the measurements for this cycle');
ncwriteatt(filename,'/flags/time','timezone',calibrationTool.timeZone);

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
ncwriteatt(filename,'/','raw_filename',calibrationTool.logFile.rawFilename);
ncwriteatt(filename,'/','raw_data_software_version',calibrationTool.logFile.raw_data_software_version);
ncwriteatt(filename,'/','calibration_version',calibrationTool.logFile.calibration_version);
ncwriteatt(filename,'/','raw_file_comment',calibrationTool.logFile.raw_file_comment);
ncwriteatt(filename,'/','comment',calibrationTool.logFile.comment);
ncwriteatt(filename,'/','filename_level1a',calibrationTool.logFile.filenameLevel1a);
ncwriteatt(filename,'/','creation_date_level1a',calibrationTool.logFile.creation_date_level1a);
ncwriteatt(filename,'/','raw_file_warning',calibrationTool.logFile.raw_file_warning);

% Global file attributes
ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(datetime('now','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','featureType','timeSeries');

ncwriteatt(filename,'/','outlier_detection',calibrationTool.outlierDectectionType);

if calibrationTool.filterByTransmittance && calibrationTool.filterByFlags
    filtering_type = 'by flags and transmittance';
elseif calibrationTool.filterByTransmittance && ~calibrationTool.filterByFlags
    filtering_type = 'by transmittance only';
elseif calibrationTool.filterByTransmittance && calibrationTool.filterByFlags
    filtering_type = 'by flags only';
else
    filtering_type = 'none';
end
ncwriteatt(filename,'/','filtering_of_calibrated_spectra',filtering_type);

if ~isempty(fieldnames(calibrationTool.labviewLog))
    if sum(isbetween([calibrationTool.labviewLog.dateTime],...
            datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day, 'TimeZone',calibrationTool.timeZone),...
            datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day+1, 'TimeZone',calibrationTool.timeZone))) > 0
        ncwriteatt(filename,'/','labview_logfile_warning','check labview log !');
    else
        ncwriteatt(filename,'/','labview_logfile_warning','clean');
    end
else
    ncwriteatt(filename,'/','labview_logfile_warning','no labview log found');
end

% Geolocation attributes
%ncwriteatt(filename,'/','data_start_date', integratedSpectra(1).first_sky_time);
%ncwriteatt(filename,'/','data_stop_date', integratedSpectra(end).last_sky_time);
if ~isnan(integratedSpectra(1).first_sky_time)
    ncwriteatt(filename,'/','data_start_date', datestr(datetime(integratedSpectra(1).first_sky_time + calibrationTool.referenceTime,'ConvertFrom','datenum','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
end
if ~isnan(integratedSpectra(end).last_sky_time)
    ncwriteatt(filename,'/','data_stop_date', datestr(datetime(integratedSpectra(end).last_sky_time + calibrationTool.referenceTime,'ConvertFrom','datenum','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
end
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
attrVal.tod = {'time of day',...
    'time_of_day',...
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

attrVal.Tb_win_corr = {'Tb_win_corr',...
    'window_corrected_brightness_temperature',...
    'K',...
    'integrated brightness temperature for this cycle, corrected for the window'};

attrVal.Tb_corr = {'Tb_corr',...
    'corrected_brightness_temperature',...
    'K',...
    'integrated brightness temperature for this cycle, corrected for window and troposphere'};

attrVal.stdTb = {'stdTb',...
    'std_brightness_temperature',...
    'K',...
    'standard deviation of brightness temperature for this cycle per channel'};

attrVal.good_channels = {'good channels',...
    'good_channels',...
    '1',...
    'Flags for the quality of the channels performed during the integration'};

attrVal.freq = {'frequency vector',...
    'frequency',...
    'Hz',...
    'frequency vector for the spectrometer'};

attrVal.if = {'intermediate frequency vector',...
    'intermediate_frequency',...
    'Hz',...
    'intermediate frequency vector for the spectrometer'};

attrVal.THot = {'THot',...
    'hot_load_temperature',...
    'K',...
    'Mean temperature of the hot load'};

% attrVal.stdTHot = {'stdTHot',...
%     'std_hot_load_temperature',...
%     'K',...
%     'standard deviation of the hot load temperature'};

attrVal.TNoise = {'noise receiver temperature',...
    'noise_temperature',...
    'K',...
    'mean noise receiver temperature'};

% attrVal.stdTNoise = {'stdTNoise',...
%     'std_noise_temperature',...
%     'K',...
%     'standard deviation of the noise receiver temperature'};
attrVal.meanTb = {'mean Tb',...
    'meanTb',...
    'K',...
    'mean brightness temperature for this cycle (without bad channel)'};

attrVal.meanStdTb = {'mean standard variation of Tb',...
    'mean_std_Tb',...
    'K',...
    'mean standard deviation of brightness temperature for this cycle (without bad channel)'};

attrVal.calibrationTime = {'calibrationTime',...
    'calibration_time',...
    'second',...
    'Time interval used for calibrating the spectra'};

attrVal.integrationTime = {'integrationTime',...
    'integration_time',...
    'second',...
    'Time interval used for integrating the spectra'};

attrVal.noiseLevel = {'std(diff(Tb))/sqrt(2)',...
    'noise_level',...
    'K',...
    'describes how noisy is the spectra'};

attrVal.meanAngleAntenna = {'mean sky angle',...
    'elevation_angle',...
    'degree',...
    'mean elevation angle of the sky observation during this cycle'};

attrVal.TRoom = {'TRoom',...
    'room_temperature',...
    'K',...
    'mean room temperature'};

attrVal.stdTRoom = {'stdTRoom',...
    'standard_room_temperature',...
    'K',...
    'standard deviation of room temperature'};

attrVal.TOut = {'TOut',...
    'outside_temperature',...
    'K',...
    'mean outside temperature'};

attrVal.TWindow = {'TWindow',...
    'window_temperature',...
    'K',...
    'mean window temperature'};

attrVal.VGunn = {'VGunn',...
    'gunn_voltage',...
    'V',...
    'mean Gunn voltage'};

% Tropospheric correction data:
attrVal.tropospheric_transmittance = {'tropospheric transmittance',...
    'tropospheric_transmittance',...
    '1',...
    'Atmospheric transmittance'};

attrVal.tropospheric_opacity = {'tropospheric opacity',...
    'tropospheric_opacity',...
    '1',...
    ''};

attrVal.tropospheric_opacity_tc = {'tropospheric opacity from tc',...
    'tropospheric_opacity_tc',...
    '1',...
    'averaged opacity derived from tipping curve measurements'};

% for Meteo data:
attrVal.air_pressure = {'air pressure',...
    'air_pressure',...
    'hPa',...
    'air pressure at the station'};

attrVal.air_temperature= {'air temperature',...
    'air_temperature',...
    'K',...
    'air temperature at the station'};

attrVal.relative_humidity = {'relative humidity',...
    'relative_humidity',...
    '1',...
    'relative humidity of the air at the station'};

attrVal.precipitation = {'precipitation',...
    'precipitation',...
    'mm',...
    'Accumulation of precipitation during the cycle (from gauge ?)'};

attrVal.numHSpectra = {'number of hot spectra',...
    'number_of_hot_spectra',...
    '1',...
    'number of hot spectra averaged together in this cycle'};

attrVal.numCSpectra = {'number of cold spectra',...
    'number_of_cold_spectra',...
    '1',...
    'number of cold spectra averaged together in this cycle'};

attrVal.numSSpectra = {'number of sky spectra',...
    'number_of_sky_spectra',...
    '1',...
    'number of sky spectra averaged together in this cycle'};

attrVal.numberOfCalSpectra = {'number of calibrated spectra',...
    'number_of_calibrated_spectra',...
    '1',...
    'number of calibrated spectra integrated during this cycle'};

for i=1:length(attrName)
    ncwriteatt(filename,'/spectrometer1/time_of_day',attrName{i},attrVal.tod{i});
    ncwriteatt(filename,'/spectrometer1/lat',attrName{i},attrVal.lat{i});
    ncwriteatt(filename,'/spectrometer1/lon',attrName{i},attrVal.lon{i});
    ncwriteatt(filename,'/spectrometer1/alt',attrName{i},attrVal.alt{i});
    ncwriteatt(filename,'/spectrometer1/azimuth_angle',attrName{i},attrVal.azimuth{i});
    %ncwriteatt(filename,'/spectrometer1/effectiveCalibrationTime',attrName{i},attrVal.effCalTime{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/Tb_win_corr',attrName{i},attrVal.Tb_win_corr{i});
    ncwriteatt(filename,'/spectrometer1/Tb_corr',attrName{i},attrVal.Tb_corr{i});
    ncwriteatt(filename,'/spectrometer1/stdTb',attrName{i},attrVal.stdTb{i});
    ncwriteatt(filename,'/spectrometer1/good_channels',attrName{i},attrVal.good_channels{i});
    ncwriteatt(filename,'/spectrometer1/frequencies',attrName{i},attrVal.freq{i});
    ncwriteatt(filename,'/spectrometer1/intermediate_freq',attrName{i},attrVal.if{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    %ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/noise_temperature',attrName{i},attrVal.TNoise{i});
    %ncwriteatt(filename,'/spectrometer1/stdTNoise',attrName{i},attrVal.stdTNoise{i});
    ncwriteatt(filename,'/spectrometer1/mean_std_Tb',attrName{i},attrVal.meanStdTb{i});
    ncwriteatt(filename,'/spectrometer1/meanTb',attrName{i},attrVal.meanTb{i});
    ncwriteatt(filename,'/spectrometer1/calibration_time',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/integration_time',attrName{i},attrVal.integrationTime{i});
    ncwriteatt(filename,'/spectrometer1/mean_sky_elevation_angle',attrName{i},attrVal.meanAngleAntenna{i});
    ncwriteatt(filename,'/spectrometer1/TRoom',attrName{i},attrVal.TRoom{i});
    ncwriteatt(filename,'/spectrometer1/stdTRoom',attrName{i},attrVal.stdTRoom{i});
    ncwriteatt(filename,'/spectrometer1/VGunn',attrName{i},attrVal.VGunn{i});
    ncwriteatt(filename,'/spectrometer1/TWindow',attrName{i},attrVal.TWindow{i});
    ncwriteatt(filename,'/spectrometer1/TOut',attrName{i},attrVal.TOut{i});
    ncwriteatt(filename,'/spectrometer1/noise_level',attrName{i},attrVal.noiseLevel{i});
    ncwriteatt(filename,'/spectrometer1/number_of_hot_spectra',attrName{i},attrVal.numHSpectra{i});
    ncwriteatt(filename,'/spectrometer1/number_of_cold_spectra',attrName{i},attrVal.numCSpectra{i});
    ncwriteatt(filename,'/spectrometer1/number_of_sky_spectra',attrName{i},attrVal.numSSpectra{i});
    ncwriteatt(filename,'/spectrometer1/number_of_calibrated_spectra',attrName{i},attrVal.numberOfCalSpectra{i});
    % Meteo
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.air_temperature{i});
    ncwriteatt(filename,'/meteo/air_pressure',attrName{i},attrVal.air_pressure{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relative_humidity{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});
    
    % Corrections
    ncwriteatt(filename,'/spectrometer1/tropospheric_transmittance',attrName{i},attrVal.tropospheric_transmittance{i});
    ncwriteatt(filename,'/spectrometer1/tropospheric_opacity',attrName{i},attrVal.tropospheric_opacity{i});
    ncwriteatt(filename,'/spectrometer1/tropospheric_opacity_tc',attrName{i},attrVal.tropospheric_opacity_tc{i});
    
end

disp(['File saved as: ' filename])
end
