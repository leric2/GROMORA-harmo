function retrievalTool = save_level1b_daily(retrievalTool,level1b)
%==========================================================================
% NAME          | save_level1a_daily.m
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
locationLevel1b=retrievalTool.level1Folder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here LOOOP in all calibration cycle and write every variable into daily vector
% of matrices

integratedSpectra=level1b.integration;

% Convert NaN's into the standard FillValue
%fields=fieldnames(integratedSpectra);
%for field = 



%initialize some variable (only the matrices)
channelId=int64(ones(length(integratedSpectra),retrievalTool.numberOfChannels)*NaN);
%Tb=ones(length(calibratedSpectra),retrievalTool.numberOfChannels)*NaN;
%frequencyVector=ones(length(calibratedSpectra),retrievalTool.numberOfChannels)*NaN;

if isfield(integratedSpectra,'errorVector')
    error=int64(ones(length(integratedSpectra),length(integratedSpectra(1).errorVector))*NaN);
    errorCalib=int64(ones(length(integratedSpectra),length(integratedSpectra(1).errorVector))*NaN);
end
for t = 1:length(integratedSpectra)
    channelId(t,:)=1:retrievalTool.numberOfChannels;
    %Tb(t,:)=calibratedSpectra(t).Tb;
    %frequencyVector(t,:)=calibratedSpectra(t).freq;
    
    if isfield(integratedSpectra,'errorVector')
        error(t,:)=1:length(integratedSpectra(1).errorVector);
        errorCalib(t,:)=integratedSpectra(t).errorVector;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now write daily level 1a file
filename=[locationLevel1b retrievalTool.instrumentName '_level1b_' retrievalTool.spectrometer '_' retrievalTool.dateStr '.nc'];
retrievalTool.filenameLevel1b=filename;

%filename=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(t) '.nc'];
%title=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(i)];

if isfile(filename)
    delete(filename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Scientific Dataset (spectrometer1,2,...)

% First create coordinates variable (enable 'netcdf4' format)
% 
nccreate(filename,'/spectrometer1/time','Dimensions',{'time',Inf},'Datatype','double','Format','netcdf4');
nccreate(filename,'/spectrometer1/channel_idx','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','int64','FillValue',-9999)

%%%%%%%%%%%%%%%%%
% Some variables to help geolocate the file:
% these are scalar variables as they do not vary in time
nccreate(filename,'/spectrometer1/lat','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/lon','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/alt','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/azimuthAngle','Dimensions',{'time',Inf},'Datatype','single','FillValue',-9999);

% some variable for better identifying the time period of the measurements
nccreate(filename,'/spectrometer1/year','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/month','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/day','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/timeOfDay','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);

nccreate(filename,'/spectrometer1/firstSkyTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/lastSkyTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);

nccreate(filename,'/spectrometer1/timeMin','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);

%%%%%%%%%%%%%%%%%
% the variables linked with the calibration    
%nccreate(filename,'/spectrometer1/effectiveCalibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/Tb','Dimensions',{'channel_idx',retrievalTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/frequencies','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/intermediateFreq','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/integrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/meanAngleAntenna','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/TRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TWindow','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_calibrated_spectra','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

% Tropospheric correction data:
nccreate(filename,'/spectrometer1/trospheric_transmittance','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/trospheric_opacity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

% Meteo Data
nccreate(filename,'/meteo/air_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/relative_humidity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/precipitation','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags dataset
if isfield(integratedSpectra,'errorVector')
    nccreate(filename,'/flags/time','Dimensions',{'time',Inf},'Datatype','double')
    nccreate(filename,'/flags/flags','Dimensions',{'flags',length(integratedSpectra(1).errorVector)},'Datatype','int64')
    
    % We input a (1xerrorVectorSize) int vector to identify the errors
    nccreate(filename,'/flags/calibration_flags','Dimensions',{'flags',length(integratedSpectra(1).errorVector),'time',Inf},'Datatype','int64','FillValue',-9999)
    %nccreate(filename,'/error/errorLabel','Dimensions',{'errorLabel',20,'time',Inf},'Datatype','char','FillValue','')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the spectrometer1 variable
% Coordinate variables, directly adding the attributes
ncwrite(filename,'/spectrometer1/time',[integratedSpectra.meanDatetime]);
ncwriteatt(filename,'/spectrometer1/time','units',retrievalTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time','calendar',retrievalTool.calendar);
ncwriteatt(filename,'/spectrometer1/time','description','mean time of THE BEGINNING of all antenna measurements for this cycle');

ncwrite(filename,'/spectrometer1/channel_idx',1:retrievalTool.numberOfChannels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some variables to help geolocate the file:
% these are scalar variables as they do not vary in time
ncwrite(filename,'/spectrometer1/lat',ones(length(integratedSpectra),1)*retrievalTool.lat);
ncwrite(filename,'/spectrometer1/lon',ones(length(integratedSpectra),1)*retrievalTool.lon);
ncwrite(filename,'/spectrometer1/alt',ones(length(integratedSpectra),1)*retrievalTool.altitude);
ncwrite(filename,'/spectrometer1/azimuthAngle',ones(length(integratedSpectra),1)*retrievalTool.azimuthAngle);

% some variable for better identifying the time period of the measurements
ncwrite(filename,'/spectrometer1/year',int64([integratedSpectra.year]));
ncwrite(filename,'/spectrometer1/month',int64([integratedSpectra.month]));
ncwrite(filename,'/spectrometer1/day',int64([integratedSpectra.day]));
ncwrite(filename,'/spectrometer1/timeOfDay',[integratedSpectra.TOD]);

ncwrite(filename,'/spectrometer1/firstSkyTime',[integratedSpectra.firstSkyTime]);  
ncwriteatt(filename,'/spectrometer1/firstSkyTime','units',retrievalTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/firstSkyTime','calendar',retrievalTool.calendar);
ncwriteatt(filename,'/spectrometer1/firstSkyTime','description','start time of the first sky measurements in this cycle');

ncwrite(filename,'/spectrometer1/lastSkyTime',[integratedSpectra.lastSkyTime]);
ncwriteatt(filename,'/spectrometer1/lastSkyTime','units',retrievalTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/lastSkyTime','calendar',retrievalTool.calendar);
ncwriteatt(filename,'/spectrometer1/lastSkyTime','description','stop time of the first sky measurements in this cycle');

ncwrite(filename,'/spectrometer1/timeMin',[integratedSpectra.timeMin]);
ncwriteatt(filename,'/spectrometer1/lastSkyTime','units',retrievalTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/lastSkyTime','calendar',retrievalTool.calendar);
ncwriteatt(filename,'/spectrometer1/lastSkyTime','description','minimum theoretical start time for this calibration cycle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main set of variables
ncwrite(filename,'/spectrometer1/Tb',vertcat(integratedSpectra.Tb)');
ncwrite(filename,'/spectrometer1/frequencies',integratedSpectra(1).freq);
ncwrite(filename,'/spectrometer1/intermediateFreq',integratedSpectra(1).if);
ncwrite(filename,'/spectrometer1/THot',[integratedSpectra.THot]);
%ncwrite(filename,'/spectrometer1/stdTHot',[integratedSpectra.stdTHot]);
ncwrite(filename,'/spectrometer1/TSys',[integratedSpectra.TSys]);
%ncwrite(filename,'/spectrometer1/stdTSys',[integratedSpectra.stdTSys]);
ncwrite(filename,'/spectrometer1/calibrationTime',[integratedSpectra.calibrationTime]);
ncwrite(filename,'/spectrometer1/integrationTime',[integratedSpectra.integrationTime]);
ncwrite(filename,'/spectrometer1/meanAngleAntenna',[integratedSpectra.meanAngleAntenna]);

ncwrite(filename,'/spectrometer1/number_calibrated_spectra',[integratedSpectra.numberOfAveragedSpectra]);

% Tropospheric correction data:
ncwrite(filename,'/spectrometer1/trospheric_transmittance',[integratedSpectra.troposphericTransmittance]);
ncwriteatt(filename,'/spectrometer1/trospheric_transmittance','method',integratedSpectra(1).transmittanceMethod);

ncwrite(filename,'/spectrometer1/trospheric_opacity',[integratedSpectra.troposphericOpacity]);
ncwriteatt(filename,'/spectrometer1/trospheric_opacity','method',integratedSpectra(1).transmittanceMethod);

% Data that are not present for every instrument
if isfield(integratedSpectra,'TempRoom')
    ncwrite(filename,'/spectrometer1/TRoom',[integratedSpectra.TempRoom]);
    ncwrite(filename,'/spectrometer1/stdTRoom',[integratedSpectra.stdTempRoom]);
else
    ncwrite(filename,'/spectrometer1/TRoom',-9999*ones(length(integratedSpectra),1));
    ncwrite(filename,'/spectrometer1/stdTRoom',-9999*ones(length(integratedSpectra),1));
end

if isfield(integratedSpectra,'TempOut')
    ncwrite(filename,'/spectrometer1/TOut',[integratedSpectra.TempOut]);
else
    ncwrite(filename,'/spectrometer1/TOut',-9999*ones(length(integratedSpectra),1));
end

if isfield(integratedSpectra,'TempWindow')
    ncwrite(filename,'/spectrometer1/TWindow',[integratedSpectra.TempWindow]);
else
    ncwrite(filename,'/spectrometer1/TWindow',-9999*ones(length(integratedSpectra),1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
ncwrite(filename,'/meteo/air_temperature',[integratedSpectra.meanAirTemperature]);
ncwrite(filename,'/meteo/relative_humidity',[integratedSpectra.meanRelativeHumidity]);
ncwrite(filename,'/meteo/precipitation',[integratedSpectra.rainAccumulation]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the flags variables
if isfield(integratedSpectra,'errorVector')
    ncwrite(filename,'/flags/time',[integratedSpectra.meanDatetime]);
    ncwriteatt(filename,'/flags/time','units',retrievalTool.meanDatetimeUnit);
    ncwriteatt(filename,'/flags/time','calendar',retrievalTool.calendar);
    ncwriteatt(filename,'/flags/time','description','mean time of the measurements for this cycle');
    
    ncwrite(filename,'/flags/flags',1:length(integratedSpectra(1).errorVector));
    ncwrite(filename,'/flags/calibration_flags',errorCalib');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originator attributes
ncwriteatt(filename,'/','title','Brightness temperature measured by ground-based radiometer');
ncwriteatt(filename,'/','location',retrievalTool.dataLocation);
ncwriteatt(filename,'/','source',retrievalTool.dataSource);
ncwriteatt(filename,'/','name',retrievalTool.PI_NAME);
ncwriteatt(filename,'/','institution',retrievalTool.PI_AFFILIATION);
ncwriteatt(filename,'/','contact',retrievalTool.PI_ADDRESS);
ncwriteatt(filename,'/','mail',retrievalTool.PI_EMAIL);

ncwriteatt(filename,'/','history','');
ncwriteatt(filename,'/','references','');
ncwriteatt(filename,'/','comment','');
%ncwriteatt(filename,'/','DATA_VARIABLES','');

ncwriteatt(filename,'/','raw_data_filename',retrievalTool.logFile.rawFilename);
ncwriteatt(filename,'/','raw_data_software_version',retrievalTool.logFile.raw_data_software_version);
ncwriteatt(filename,'/','calibrated_version',retrievalTool.logFile.calibration_version);
ncwriteatt(filename,'/','raw_file_comment',retrievalTool.logFile.raw_file_comment);
ncwriteatt(filename,'/','comment',retrievalTool.logFile.comment);
ncwriteatt(filename,'/','filename_level1a',retrievalTool.logFile.filenameLevel1a);
ncwriteatt(filename,'/','creation_date_level1a',retrievalTool.logFile.creation_date_level1a);
ncwriteatt(filename,'/','raw_file_warning',retrievalTool.logFile.raw_file_warning);

% Global file attributes
ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(now,'yyyymmddTHHMMSSZ'));   %TODO
ncwriteatt(filename,'/','featureType','timeSeries');

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
    
    %ncwriteatt(filename,'/flags/calibration_flags','flag_meanings','sufficientNumberOfIndices systemTemperatureOK hotAngleRemoved coldAngleRemoved antennaAngleRemoved LN2SensorsOK LN2LevelOK hotLoadOK FFT_adc_overload_OK');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_1','sufficientNumberOfIndices');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_2','systemTemperatureOK');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_3','hotAngleRemoved');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_4','coldAngleRemoved');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_5','antennaAngleRemoved');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_6','LN2SensorsOK');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_7','LN2LevelOK');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_8','hotLoadOK');
    ncwriteatt(filename,'/flags/calibration_flags','errorCode_9','FFT_adc_overload_OK');
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
    'calibrated brightness temperature for this cycle'};

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
attrVal.trospheric_transmittance = {'tropospheric transmittance',...
    '',...
    '1',...
    ''};

attrVal.trospheric_opacity = {'tropospheric opacity',...
    '',...
    '1',...
    ''};

% Meteo Data
% Tropospheric correction data:
attrVal.air_temperature = {'air temperature',...
    'air_temperature',...
    'K',...
    ''};

attrVal.relative_humidity = {'relative humidity',...
    'relative_humidity',...
    '1',...
    ''};

attrVal.precipitation = {'rain accumulation',...
    '',...
    'mm',...
    ''};

for i=1:length(attrName)
    ncwriteatt(filename,'/spectrometer1/timeOfDay',attrName{i},attrVal.tod{i});
    ncwriteatt(filename,'/spectrometer1/lat',attrName{i},attrVal.lat{i});
    ncwriteatt(filename,'/spectrometer1/lon',attrName{i},attrVal.lon{i});
    ncwriteatt(filename,'/spectrometer1/alt',attrName{i},attrVal.alt{i});
    ncwriteatt(filename,'/spectrometer1/azimuthAngle',attrName{i},attrVal.azimuth{i});
    %ncwriteatt(filename,'/spectrometer1/effectiveCalibrationTime',attrName{i},attrVal.effCalTime{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/frequencies',attrName{i},attrVal.freq{i});
    ncwriteatt(filename,'/spectrometer1/intermediateFreq',attrName{i},attrVal.if{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    %ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/TSys',attrName{i},attrVal.TSys{i});
    %ncwriteatt(filename,'/spectrometer1/stdTSys',attrName{i},attrVal.stdTSys{i});
    ncwriteatt(filename,'/spectrometer1/calibrationTime',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/integrationTime',attrName{i},attrVal.integrationTime{i});
    ncwriteatt(filename,'/spectrometer1/meanAngleAntenna',attrName{i},attrVal.meanAngleAntenna{i});
    ncwriteatt(filename,'/spectrometer1/TRoom',attrName{i},attrVal.TRoom{i});
    ncwriteatt(filename,'/spectrometer1/stdTRoom',attrName{i},attrVal.stdTRoom{i});
    ncwriteatt(filename,'/spectrometer1/TWindow',attrName{i},attrVal.TWindow{i});
    ncwriteatt(filename,'/spectrometer1/TOut',attrName{i},attrVal.TOut{i});
    
    % Meteo
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.air_temperature{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relative_humidity{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});
    
    
    % Corrections
    ncwriteatt(filename,'/spectrometer1/trospheric_transmittance',attrName{i},attrVal.trospheric_transmittance{i});
    ncwriteatt(filename,'/spectrometer1/trospheric_opacity',attrName{i},attrVal.trospheric_opacity{i});
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding debug (optionnal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if retrievalTool.saveAllCycles
    % initialize matrices
    nClean=0;
    for t = 1:length(integratedSpectra)
        nClean=nClean+integratedSpectra(t).numberOfCleanAntennaAngle;
    end
    
    channel_idx=int64(ones(nClean,retrievalTool.numberOfChannels)*NaN);
    TbAll=single(ones(nClean,retrievalTool.numberOfChannels)*NaN);
    cycleNumber=int64(ones(nClean,1)*NaN);
    cycleId=int64(1:nClean);
    
    counter=1;
    for t = 1:length(integratedSpectra)
        for i =1:integratedSpectra(t).numberOfCleanAntennaAngle
            TbAll(counter,:)=integratedSpectra(t).TbAll(i,:);
            channel_idx(counter,:)=1:retrievalTool.numberOfChannels;
            cycleNumber(counter)=t;
            counter=counter+1;
        end
    end
    
    % Group for debugging variables:
    nccreate(filename,'/debug/cycle_id','Dimensions',{'cycle_id',Inf},'Datatype','int64')
    nccreate(filename,'/debug/channel_idx','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','int64','FillValue',-9999)
    
    nccreate(filename,'/debug/calibration_cycle_number','Dimensions',{'cycle_id',Inf},'Datatype','int64','FillValue',-9999)
    nccreate(filename,'/debug/Tb_all','Dimensions',{'channel_idx',retrievalTool.numberOfChannels,'cycle_id',Inf},'Datatype','single','FillValue',-9999)
    
    % Writing the debug variables
    ncwrite(filename,'/debug/cycle_id',cycleId);
    
    ncwrite(filename,'/debug/channel_idx',1:retrievalTool.numberOfChannels);
    ncwrite(filename,'/debug/calibration_cycle_number',cycleNumber');
    ncwrite(filename,'/debug/Tb_all',TbAll');
end

disp(['File saved as: ' filename])
end
