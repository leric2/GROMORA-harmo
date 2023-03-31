function calibrationTool = save_l1a_LN2(calibrationTool,logFile,calibratedSpectra,warningLevel0)
%==========================================================================
% NAME      | save_level1a_daily.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 01.2020
%           |
% ABSTRACT  | Function saving the calibrated spectra for a complete day
%           | in the form of netCDF4 file. The file contains one entry
%           | per calibration cylce (typically 10 minutes) on a time 
%           | coordinate as well as all necessery metadata to identify
%           | and re-use the data.
%           |
%           |
% ARGUMENTS | INPUTS: 1. calibrationTool:
%           |           - calType
%           |           - level1Folder
%           |           - instrumentName
%           |           - spectrometer
%           |           - dateStr
%           |           - numberOfChannels
%           |           - meanDatetimeUnit
%           |           - calendar
%           |           - lat, lon, altitude, azimuthAngle
%           |           - dataLocation, dataSource
%           |           - PI_NAME, PI_AFFILIATION, PI_ADDRESS, PI_EMAIL
%           |           - numberOfSpectrometer
%           |           - labviewLog (optional)
%           |         2. logFile
%           |         3. calibratedSpectra
%           |         4. warningLevel0
%           |         
%           | OUTPUTS: 1. calibrationTool with added field
%           |             "filenameLevel1a"
%           |
% SAVE      | Level1a netCDF4 containing the following groups:
%           |   - '/' for global attributes
%           |   - '/spectrometer1/' for spectrometer dataset
%           |   - '/flags/'
%           |   - '/meteo/'
%           |  
%==========================================================================
% Filename and location for DAILY netCDF file
disp('starting save file')

if strcmp(calibrationTool.calType,'debug')
    filename=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '_debug.nc'];
else
    filename=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr calibrationTool.extraName '.nc'];
end
calibrationTool.filenameLevel1a=filename;

% Rewrite if already existing
if  exist(filename,'file')
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

for i = 1:calibrationTool.spectrometerQuantity
    
    N_spectro = sprint('/spectrometer%d', i);
    nccreate(filename,N_spectro+'/lat','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
    nccreate(filename,N_spectro+'/lon','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
    nccreate(filename,N_spectro+'/alt','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
    nccreate(filename,N_spectro+'/azimuth_angle','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
    
    % time variables
    nccreate(filename,N_spectro+'/MJD2K','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,N_spectro+'/year','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
    nccreate(filename,N_spectro+'/month','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
    nccreate(filename,N_spectro+'/day','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
    nccreate(filename,N_spectro+'/time_of_day','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,N_spectro+'/first_sky_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,N_spectro+'/last_sky_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,N_spectro+'/time_min','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
end
%%%%%%%%%%%%%%%%%
% Calibration variables   
%nccreate(filename,'/spectrometer1/effectiveCalibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
if calibrationTool.spectrometerQuantity == 1
    nccreate(filename,'/spectrometer1/Tb_line','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/spectrometer1/stdTb_line','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/spectrometer1/frequency','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/spectrometer1/intermediate_freq','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    
    nccreate(filename,'/spectrometer1/Tb_ref','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/spectrometer1/stdTb_ref','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)

    nccreate(filename,'/spectrometer1/Tb_cold','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/spectrometer1/stdTb_cold','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)

    nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

    nccreate(filename,'/spectrometer1/TSys','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/spectrometer1/stdTSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/spectrometer1/calibration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

elseif calibrationTool.spectrometerQuantity > 1
    for i = 1:calibrationTool.spectrometerQuantity
      N_spectro = sprint('spectrometer%d', i);
      if i == 1
          chanNos = calibrationTool.spectrometer1.numberOfChannels;
      elseif i == 2
          chanNos = calibrationTool.spectrometer2.numberOfChannels;
      end

      nccreate(filename,N_spectro+'/Tb','Dimensions',{'channel_idx',chanNos,'time',Inf},'Datatype','double','FillValue',-9999);
      nccreate(filename,N_spectro+'/stdTb','Dimensions',{'channel_idx',chanNos,'time',Inf},'Datatype','double','FillValue',-9999);
      nccreate(filename,N_spectro+'/frequency','Dimensions',{'channel_idx',chanNos},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/intermediate_freq','Dimensions',{'channel_idx',chanNos},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/meanTb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/mean_std_Tb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/Tsys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/stdTSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
      nccreate(filename,N_spectro+'/calibration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scientific Dataset (spectrometer1,2,...)

% Coordinate variables, directly with their attributes
if calibrationTool.spectrometerQuantity == 1
    ncwrite(filename,'/spectrometer1/channel_idx',1:calibrationTool.numberOfChannels);
    ncwriteatt(filename,'/spectrometer1/channel_idx','description','index of the spectrometer channels, from 1 to N (number of channels)');
    ncwrite(filename,'/spectrometer1/lat',ones(length(calibratedSpectra),1)*calibrationTool.lat);
    ncwrite(filename,'/spectrometer1/lon',ones(length(calibratedSpectra),1)*calibrationTool.lon);
    ncwrite(filename,'/spectrometer1/alt',ones(length(calibratedSpectra),1)*calibrationTool.altitude);
    ncwrite(filename,'/spectrometer1/azimuth_angle',ones(length(calibratedSpectra),1)*calibrationTool.azimuthAngle);
    
    % Time variables
    ncwrite(filename,'/spectrometer1/year',int64([calibratedSpectra.year]));
    ncwriteatt(filename,'/spectrometer1/year','description','year of the measurement as integer');
    
    ncwrite(filename,'/spectrometer1/month',int64([calibratedSpectra.month]));
    ncwriteatt(filename,'/spectrometer1/month','description','month of the measurement as integer');
    
    ncwrite(filename,'/spectrometer1/day',int64([calibratedSpectra.day]));
    ncwriteatt(filename,'/spectrometer1/day','description','day of the month as integer');
    
    ncwrite(filename,'/spectrometer1/time_of_day',[calibratedSpectra.timeOfDay]);
    
    ncwrite(filename,'/spectrometer1/first_sky_time',[calibratedSpectra.firstSkyTime]);  
    ncwriteatt(filename,'/spectrometer1/first_sky_time','units',calibrationTool.meanDatetimeUnit);
    ncwriteatt(filename,'/spectrometer1/first_sky_time','calendar',calibrationTool.calendar);
    ncwriteatt(filename,'/spectrometer1/first_sky_time','description','time of the first sky measurements in this calibration cycle');
    
    ncwrite(filename,'/spectrometer1/last_sky_time',[calibratedSpectra.lastSkyTime]);
    ncwriteatt(filename,'/spectrometer1/last_sky_time','units',calibrationTool.meanDatetimeUnit);
    ncwriteatt(filename,'/spectrometer1/last_sky_time','calendar',calibrationTool.calendar);
    ncwriteatt(filename,'/spectrometer1/last_sky_time','description','time of the last sky measurements in this calibration cycle');
    
    ncwrite(filename,'/spectrometer1/time_min',[calibratedSpectra.timeMin]);
    ncwriteatt(filename,'/spectrometer1/time_min','units',calibrationTool.meanDatetimeUnit);
    ncwriteatt(filename,'/spectrometer1/time_min','calendar',calibrationTool.calendar);
    ncwriteatt(filename,'/spectrometer1/time_min','description','minimum theoretical start time for this calibration cycle');

    ncwrite(filename,'/spectrometer1/Tb_line',vertcat(calibratedSpectra.Tb_line)');
    ncwrite(filename,'/spectrometer1/Tb_ref',vertcat(calibratedSpectra.Tb_ref)');
    ncwrite(filename,'/spectrometer1/Tb_cold',vertcat(calibratedSpectra.TCold)');

    ncwrite(filename,'/spectrometer1/THot',vertcat(calibratedSpectra.THot)');
    ncwrite(filename,'/spectrometer1/TSys',vertcat(calibratedSpectra.TSys)');
    ncwrite(filename,'/spectrometer1/frequency',calibratedSpectra(1).freq);

    ncwrite(filename,'/spectrometer1/meanAngleLine',vertcat(calibratedSpectra.meanAngleLine)');
    ncwrite(filename,'/spectrometer1/meanAngleRef',vertcat(calibratedSpectra.meanAngleRef)');


elseif calibrationTool.spectrometerQuantity > 1
    for i = 1:calibrationTool.spectrometerQuantity
      N_spectro = sprint('spectrometer%d', i);
      if i == 1
          stChan = 1;
          chanNos = calibrationTool.spectrometer1.numberOfChannels;
      elseif i == 2
          stchan = calibrationTool.spectrometer1.numberOfChannels+1;
          chanNos = calibrationTool.spectrometer2.numberOfChannels;
      end
      ncwrite(filename,N_spectro+'/channel_idx',1:calibrationTool.numberOfChannels);
      ncwriteatt(filename,N_spectro+'/channel_idx','description','index of the spectrometer channels, from 1 to N (number of channels)');
      ncwrite(filename,N_spectro+'/lat',ones(length(calibratedSpectra),1)*calibrationTool.lat);
      ncwrite(filename,N_spectro+'/lon',ones(length(calibratedSpectra),1)*calibrationTool.lon);
      ncwrite(filename,N_spectro+'/alt',ones(length(calibratedSpectra),1)*calibrationTool.altitude);
      ncwrite(filename,N_spectro+'/azimuth_angle',ones(length(calibratedSpectra),1)*calibrationTool.azimuthAngle);
      
      % Time variables
      ncwrite(filename,N_spectro+'/year',int64([calibratedSpectra.year]));
      ncwriteatt(filename,N_spectro+'/year','description','year of the measurement as integer');
        
      ncwrite(filename,N_spectro+'/month',int64([calibratedSpectra.month]));
      ncwriteatt(filename,N_spectro+'/month','description','month of the measurement as integer');
        
      ncwrite(filename,N_spectro+'/day',int64([calibratedSpectra.day]));
      ncwriteatt(filename,N_spectro+'/day','description','day of the month as integer');
       
      ncwrite(filename,N_spectro+'/time_of_day',[calibratedSpectra.timeOfDay]);
      
      ncwrite(filename,N_spectro+'/first_sky_time',[calibratedSpectra.firstSkyTime]);  
      ncwriteatt(filename,N_spectro+'/first_sky_time','units',calibrationTool.meanDatetimeUnit);
      ncwriteatt(filename,N_spectro+'/first_sky_time','calendar',calibrationTool.calendar);
      ncwriteatt(filename,N_spectro+'/first_sky_time','description','time of the first sky measurements in this calibration cycle');
        
      ncwrite(filename,N_spectro+'/last_sky_time',[calibratedSpectra.lastSkyTime]);
      ncwriteatt(filename,N_spectro+'/last_sky_time','units',calibrationTool.meanDatetimeUnit);
      ncwriteatt(filename,N_spectro+'/last_sky_time','calendar',calibrationTool.calendar);
      ncwriteatt(filename,N_spectro+'/last_sky_time','description','time of the last sky measurements in this calibration cycle');
        
      ncwrite(filename,N_spectro+'/time_min',[calibratedSpectra.timeMin]);
      ncwriteatt(filename,N_spectro+'/time_min','units',calibrationTool.meanDatetimeUnit);
      ncwriteatt(filename,N_spectro+'time_min','calendar',calibrationTool.calendar);
      ncwriteatt(filename,N_spectro+'/time_min','description','minimum theoretical start time for this calibration cycle');

      ncwrite(filename,N_spectro+'/Tb_line',vertcat(calibratedSpectra.Tb_line(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/Tb_ref', vertcat(calibratedSpectra.Tb_ref(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/Tb_cold',vertcat(calibratedSpectra.TCold(stChan:chanNos,:))');

      ncwrite(filename,N_spectro+'/Tb_line',vertcat(calibratedSpectra.Tb_line(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/Tb_ref',vertcat(calibratedSpectra.Tb_ref(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/Tb_cold',vertcat(calibratedSpectra.TCold(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/THot',vertcat(calibratedSpectra.THot(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/TSys',vertcat(calibratedSpectra.TSys(stChan:chanNos,:))');

      ncwrite(filename,N_spectro+'/stdTSys',vertcat(calibratedSpectra.stdTSys(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/frequency',calibratedSpectra(1).freq(stChan:chanNos,:)');

      ncwrite(filename,N_spectro+'/meanAngleLine',vertcat(calibratedSpectra.meanAngleLine(stChan:chanNos,:))');
      ncwrite(filename,N_spectro+'/meanAngleRef',calibratedSpectra(1).meanAngleRef(stChan:chanNos,:)');

    end
end

%%%%%%%%%%%%%%%%%
% Calibration variables

%ncwrite(filename,'/spectrometer1/stdTb',vertcat(calibratedSpectra.stdTb)');
if isfield(calibratedSpectra,'stdTb')
    ncwrite(filename,'/spectrometer1/stdTb',vertcat(calibratedSpectra.stdTb)');
else
    ncwrite(filename,'/spectrometer1/stdTb',-9999*ones(length(calibratedSpectra),1));
end

%ncwrite(filename,'/spectrometer1/mean_std_Tb',[calibratedSpectra.meanStdTb]);

if isfield(calibratedSpectra,'meanStdTb')
    ncwrite(filename,'/spectrometer1/mean_std_Tb',[calibratedSpectra.meanStdTb]);
else
    ncwrite(filename,'/spectrometer1/mean_std_Tb',-9999*ones(length(calibratedSpectra),1));
end

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
ncwriteatt(filename,'/','raw_filename',logFile.file);
if isfield(logFile,'SW_version')
    ncwriteatt(filename,'/','raw_data_software_version',num2str(logFile.SW_version(1)));
else
    ncwriteatt(filename,'/','raw_data_software_version','unknown');
end
ncwriteatt(filename,'/','calibration_version',calibratedSpectra(1).calibrationVersion);
ncwriteatt(filename,'/','raw_file_comment',logFile.comment);

ncwriteatt(filename,'/','raw_file_warning',warningLevel0);

ncwriteatt(filename,'/','outlier_detection',calibrationTool.outlierDectectionType);

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
%ncwriteatt(filename,'/','data_start_date',calibratedSpectra(1).firstSkyTime);
%ncwriteatt(filename,'/','data_stop_date',calibratedSpectra(end).lastSkyTime);
if ~isnan(calibratedSpectra(1).firstSkyTime)
    ncwriteatt(filename,'/','data_start_date', datestr(datetime(calibratedSpectra(1).firstSkyTime + calibrationTool.referenceTime,'ConvertFrom','datenum','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
end
if ~isnan(calibratedSpectra(end).lastSkyTime)
    ncwriteatt(filename,'/','data_stop_date', datestr(datetime(calibratedSpectra(end).lastSkyTime + calibrationTool.referenceTime,'ConvertFrom','datenum','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
end

ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(datetime('now','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','featureType','timeSeries');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attributes for groups and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the attributes of flags group

if isfield(calibratedSpectra,'errorVector')
    ncwriteatt(filename,'/flags','description','Each spectra is associated with a flag vector. The order and meaning of the flags are described in its attributes');
    ncwriteatt(filename,'/flags','number_of_flags',lenErrorVect);

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

attrVal.tod = {'Time of day',...
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


attrVal.Tb = {'Tb',...
    'brightness_temperature',...
    'K',...
    'calibrated brightness temperature for this cycle'};

attrVal.stdTb = {'standard variation of Tb',...
    'std_Tb',...
    'K',...
    'standard deviation of brightness temperature for this cycle per channel'};

attrVal.freq = {'frequency vector',...
    'frequency',...
    'Hz',...
    'frequency vector for the spectrometer'};

attrVal.meanStdTb = {'mean standard variation of Tb',...
    'mean_std_Tb',...
    'K',...
    'mean standard deviation of brightness temperature for this cycle (without bad channel)'};

attrVal.meanTb = {'mean Tb',...
    'meanTb',...
    'K',...
    'mean brightness temperature for this cycle (without bad channel)'};

attrVal.if = {'intermediate frequency vector',...
    'intermediate_frequency',...
    'Hz',...
    'intermediate frequency vector for the spectrometer'};

attrVal.THot = {'THot',...
    'hot_load_temperature',...
    'K',...
    'Mean temperature of the hot load'};

attrVal.stdTHot = {'stdTHot',...
    'std_hot_load_temperature',...
    'K',...
    'standard deviation of the hot load temperature'};

attrVal.TNoise = {'noise receiver temperature',...
    'noise_temperature',...
    'K',...
    'mean noise receiver temperature'};

attrVal.stdTNoise = {'standard deviation of noise receiver temperature',...
    'std_noise_temperature',...
    'K',...
    'standard deviation of the noise receiver temperature'};

attrVal.calibrationTime = {'calibrationTime',...
    'calibration_time',...
    'second',...
    'Time interval used for calibrating the spectra'};

attrVal.noiseLevel = {'std(diff(Tb))/sqrt(2)',...
    'noise_level',...
    'K',...
    'describes how noisy is the spectra'};

attrVal.meanAngleAntenna = {'mean sky angle',...
    'elevation_angle',...
    'degree',...
    'mean elevation angle of the sky observation during this cycle'};

attrVal.meanOpacity = {'meanOpaciy',...
    '',...
    '',...
    'mean zenith opacity of the atmosphere during this cycle'};

attrVal.sigma = {'sigma',...
    '',...
    'K',...
    'noise temperature of the integrated spectrum'};

attrVal.A = {'A',...
    '',...
    '',...
    'the equivalent transmission of the reference absorber'};

attrVal.a = {'a',...
    '',...
    '',...
    'the correction coefficient for the troposphere and ref absorber'};

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

attrVal.stdVGunn = {'stdVGunn',...
    'standard_gunn_voltage',...
    'V',...
    'standard deviation of Gunn voltage'};

% for Meteo data:
attrVal.airP = {'air pressure',...
    'air_pressure',...
    'hPa',...
    'air pressure at the station'};

attrVal.airT = {'air temperature',...
    'air_temperature',...
    'K',...
    'air temperature at the station'};

attrVal.relH = {'relative humidity',...
    'relative_humidity',...
    '%',...
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

attrVal.meanColdCounts = {'mean FFTS cold counts',...
    'mean_cold_counts',...
    '1',...
    'mean raw FFTS counts on cold load during this cycle'};

attrVal.meanHotCounts = {'mean FFTS hot counts',...
    'mean_hot_counts',...
    '1',...
    'mean raw FFTS counts on hot load during this cycle'};

attrVal.meanSkyCounts = {'mean FFTS sky counts',...
    'mean_sky_counts',...
    '1',...
    'mean raw FFTS counts for sky observation during this cycle'};

attrVal.Tb_line = {'Mean brightness temperature sky observation',...
    'Tb_line',...
    'K',...
    'mean calibrated brightness temperature for sky observation during this cycle '};

attrVal.Tb_ref = {'Mean brightness temperature reference observation',...
    'Tb_ref',...
    'K',...
    'mean calibrated brightness temperature for reference observation during this cycle '};

attrVal.Tb_cold = {'Mean brightness temperature cold sky observation',...
    'Tb_cold',...
    'K',...
    'mean calibrated brightness temperature for cold sky observation during this cycle '};

attrVal.Tb_cold = {'Mean brightness temperature cold sky observation',...
    'Tb_cold',...
    'K',...
    'mean calibrated brightness temperature for cold sky observation during this cycle '};

attrVal.THot = {'Mean temperature of the ambient absorber',...
    'THot',...
    'K',...
    'mean physical temperature of the hot reference target during this cycle '};

attrVal.TSys = {'mean system noise temperature',...
    'TSys',...
    'K',...
    ['Mean system noise temperature at each frequency throughout ' ...
    'the cycle ']};

attrVal.stdTSys = {'Standard deviation of the system noise temperature',...
    'stdTSys',...
    'K',...
    ['Standard deviation of the system noise temperature at each frequency throughout ' ...
    'the cycle ']};

attrVal.meanAngleLine = {'Mean elevation angle of the line observation',...
    'meanAngleLine',...
    'degrees',...
    'Mean elevation angle of the sky observation across the cycle '};

attrVal.meanAngleRef = {'Mean elevation angle of the reference observation',...
    'meanAngleRef',...
    'degrees',...
    'Mean elevation angle of the reference observation across the cycle'};

% Ugly and open to suggestion
for i=1:length(attrName)
  if calibrationTool.spectrometerQuantity == 1

    ncwriteatt(filename,'/spectrometer1/time_of_day',attrName{i},attrVal.tod{i});
    ncwriteatt(filename,'/spectrometer1/lat',attrName{i},attrVal.lat{i});
    ncwriteatt(filename,'/spectrometer1/lon',attrName{i},attrVal.lon{i});
    ncwriteatt(filename,'/spectrometer1/alt',attrName{i},attrVal.alt{i});
    ncwriteatt(filename,'/spectrometer1/azimuth_angle',attrName{i},attrVal.azimuth{i});

    %ncwriteatt(filename,'/spectrometer1/effectiveCalibrationTime',attrName{i},attrVal.effCalTime{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/stdTb',attrName{i},attrVal.stdTb{i});
    ncwriteatt(filename,'/spectrometer1/mean_std_Tb',attrName{i},attrVal.meanStdTb{i});
    ncwriteatt(filename,'/spectrometer1/meanTb',attrName{i},attrVal.meanTb{i});
    ncwriteatt(filename,'/spectrometer1/frequency',attrName{i},attrVal.freq{i});
    ncwriteatt(filename,'/spectrometer1/intermediate_freq',attrName{i},attrVal.if{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/noise_temperature',attrName{i},attrVal.TNoise{i});
    ncwriteatt(filename,'/spectrometer1/std_dev_noise_temperature',attrName{i},attrVal.stdTNoise{i});
    ncwriteatt(filename,'/spectrometer1/calibration_time',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/mean_opacity',attrName{i},attrVal.meanOpacity{i});
    ncwriteatt(filename,'/spectrometer1/sigma',attrName{i},attrVal.sigma{i});
    ncwriteatt(filename,'/spectrometer1/A',attrName{i},attrVal.A{i});
    ncwriteatt(filename,'/spectrometer1/a',attrName{i},attrVal.a{i});
    ncwriteatt(filename,'/spectrometer1/TRoom',attrName{i},attrVal.TRoom{i});
    ncwriteatt(filename,'/spectrometer1/stdTRoom',attrName{i},attrVal.stdTRoom{i});
    ncwriteatt(filename,'/spectrometer1/VGunn',attrName{i},attrVal.VGunn{i});
    ncwriteatt(filename,'/spectrometer1/stdVGunn',attrName{i},attrVal.stdVGunn{i});
    ncwriteatt(filename,'/spectrometer1/TWindow',attrName{i},attrVal.TWindow{i});
    ncwriteatt(filename,'/spectrometer1/TOut',attrName{i},attrVal.TOut{i});
    ncwriteatt(filename,'/spectrometer1/noise_level',attrName{i},attrVal.noiseLevel{i});
    ncwriteatt(filename,'/spectrometer1/number_of_hot_spectra',attrName{i},attrVal.numHSpectra{i});
    ncwriteatt(filename,'/spectrometer1/number_of_cold_spectra',attrName{i},attrVal.numCSpectra{i});
    ncwriteatt(filename,'/spectrometer1/number_of_sky_spectra',attrName{i},attrVal.numSSpectra{i});
    ncwriteatt(filename,'/spectrometer1/mean_cold_counts',attrName{i},attrVal.meanColdCounts{i});
    ncwriteatt(filename,'/spectrometer1/mean_hot_counts',attrName{i},attrVal.meanHotCounts{i});
    ncwriteatt(filename,'/spectrometer1/mean_sky_counts',attrName{i},attrVal.meanSkyCounts{i});

    ncwriteatt(filename,'/spectrometer1/Tb_line',attrName{i},attrVal.Tb_line{i});
    ncwriteatt(filename,'/spectrometer1/Tb_ref',attrName{i},attrVal.Tb_ref{i});
    ncwriteatt(filename,'/spectrometer1/Tb_cold',attrName{i},attrVal.Tb_cold{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    ncwriteatt(filename,'/spectrometer1/stdTSys',attrName{i},attrVal.stdTSys{i});  
    ncwriteatt(filename,'/spectrometer1/meanAngleRef',attrName{i},attrVal.meanAngleRef{i});
    ncwriteatt(filename,'/spectrometer1/meanAngleLine',attrName{i},attrVal.meanAngleLine{i});
    
    % Meteo attr
    ncwriteatt(filename,'/meteo/air_pressure',attrName{i},attrVal.airP{i});
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.airT{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relH{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});

  elseif calibrationTool.spectrometerQuantity > 1
    for i = 1:calibrationTool.spectrometerQuantity
        N_spectro = sprint('spectrometer%d', i);
        if i == 1
            chanNos = calibrationTool.spectrometer1.numberOfChannels;
        elseif i == 2
            chanNos = calibrationTool.spectrometer2.numberOfChannels;
        end

        ncwriteatt(filename,N_spectro+'/time_of_day',attrName{i},attrVal.tod{i});
        ncwriteatt(filename,N_spectro+'/lat',attrName{i},attrVal.lat{i});
        ncwriteatt(filename,N_spectro+'/lon',attrName{i},attrVal.lon{i});
        ncwriteatt(filename,N_spectro+'/alt',attrName{i},attrVal.alt{i});
        ncwriteatt(filename,N_spectro+'/azimuth_angle',attrName{i},attrVal.azimuth{i});
        
        %ncwriteatt(filename,'/spectrometer1/effectiveCalibrationTime',attrName{i},attrVal.effCalTime{i});
        ncwriteatt(filename,N_spectro+'/Tb',attrName{i},attrVal.Tb{i});
        ncwriteatt(filename,N_spectro+'/stdTb',attrName{i},attrVal.stdTb{i});
        ncwriteatt(filename,N_spectro+'/mean_std_Tb',attrName{i},attrVal.meanStdTb{i});
        ncwriteatt(filename,N_spectro+'/meanTb',attrName{i},attrVal.meanTb{i});
        ncwriteatt(filename,N_spectro+'/frequency',attrName{i},attrVal.freq{i});
        ncwriteatt(filename,N_spectro+'/intermediate_freq',attrName{i},attrVal.if{i});
        ncwriteatt(filename,N_spectro+'/THot',attrName{i},attrVal.THot{i});
        ncwriteatt(filename,N_spectro+'/stdTHot',attrName{i},attrVal.stdTHot{i});
        ncwriteatt(filename,N_spectro+'/noise_temperature',attrName{i},attrVal.TNoise{i});
        ncwriteatt(filename,N_spectro+'/std_dev_noise_temperature',attrName{i},attrVal.stdTNoise{i});
        ncwriteatt(filename,N_spectro+'/calibration_time',attrName{i},attrVal.calibrationTime{i});
        ncwriteatt(filename,N_spectro+'/mean_opacity',attrName{i},attrVal.meanOpacity{i});
        ncwriteatt(filename,N_spectro+'/sigma',attrName{i},attrVal.sigma{i});
        ncwriteatt(filename,N_spectro+'/A',attrName{i},attrVal.A{i});
        ncwriteatt(filename,N_spectro+'/a',attrName{i},attrVal.a{i});
        ncwriteatt(filename,N_spectro+'/TRoom',attrName{i},attrVal.TRoom{i});
        ncwriteatt(filename,N_spectro+'/stdTRoom',attrName{i},attrVal.stdTRoom{i});
        ncwriteatt(filename,N_spectro+'/VGunn',attrName{i},attrVal.VGunn{i});
        ncwriteatt(filename,N_spectro+'/stdVGunn',attrName{i},attrVal.stdVGunn{i});
        ncwriteatt(filename,N_spectro+'/TWindow',attrName{i},attrVal.TWindow{i});
        ncwriteatt(filename,N_spectro+'/TOut',attrName{i},attrVal.TOut{i});
        ncwriteatt(filename,N_spectro+'/noise_level',attrName{i},attrVal.noiseLevel{i});
        ncwriteatt(filename,N_spectro+'/number_of_hot_spectra',attrName{i},attrVal.numHSpectra{i});
        ncwriteatt(filename,N_spectro+'/number_of_cold_spectra',attrName{i},attrVal.numCSpectra{i});
        ncwriteatt(filename,N_spectro+'/number_of_sky_spectra',attrName{i},attrVal.numSSpectra{i});
        ncwriteatt(filename,N_spectro+'/mean_cold_counts',attrName{i},attrVal.meanColdCounts{i});
        ncwriteatt(filename,N_spectro+'/mean_hot_counts',attrName{i},attrVal.meanHotCounts{i});
        ncwriteatt(filename,N_spectro+'/mean_sky_counts',attrName{i},attrVal.meanSkyCounts{i});
        
        ncwriteatt(filename,N_spectro+'/Tb_line',attrName{i},attrVal.Tb_line{i});
        ncwriteatt(filename,N_spectro+'/Tb_ref',attrName{i},attrVal.Tb_ref{i});
        ncwriteatt(filename,N_spectro+'/Tb_cold',attrName{i},attrVal.Tb_cold{i});
        ncwriteatt(filename,N_spectro+'/THot',attrName{i},attrVal.THot{i});
        ncwriteatt(filename,N_spectro+'/stdTSys',attrName{i},attrVal.stdTSys{i});  
        ncwriteatt(filename,N_spectro+'/meanAngleRef',attrName{i},attrVal.meanAngleRef{i});
        ncwriteatt(filename,N_spectro+'/meanAngleLine',attrName{i},attrVal.meanAngleLine{i});
        
    % Meteo attr
    ncwriteatt(filename,'/meteo/air_pressure',attrName{i},attrVal.airP{i});
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.airT{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relH{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});
 
end



disp(['File saved as: ' filename])
end
