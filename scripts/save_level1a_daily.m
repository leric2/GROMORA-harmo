function feedback = save_level1a_daily(retrievalTool,log,calibratedSpectra)
%==========================================================================
% NAME          | 
% TYPE          |
% AUTHOR(S)     |
% CREATION      |
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

filename=[locationLevel1a retrievalTool.instrumentName '_level1a_' retrievalTool.dateStr '.nc'];
%filename=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(t) '.nc'];
%title=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(i)];

if exist(filename)
    delete(filename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the different dataset
% Scientific Dataset (spectrometer1)
% First create coordinates variable (enable 'netcdf4' format)
nccreate(filename,'/spectrometer1/time','Dimensions',{'time',Inf},'Datatype','double','Format','netcdf4');
nccreate(filename,'/spectrometer1/channel_idx','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','int64','FillValue',-9999)

% 
nccreate(filename,'/spectrometer1/year','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/month','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/day','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999);
nccreate(filename,'/spectrometer1/timeOfDay','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);

% Scalar variables
nccreate(filename,'/spectrometer1/lat','Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/lon','Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/alt','Datatype','single','FillValue',-9999);
nccreate(filename,'/spectrometer1/azimuthAngle','Datatype','single','FillValue',-9999);

%nccreate(filename,'/spectrometer1/startTime','Dimensions',{'time',Inf},'Datatype','float','FillValue','NA');
%nccreate(filename,'/spectrometer1/stopTime','Dimensions',{'time',Inf},'Datatype','char','FillValue','NA');

nccreate(filename,'/spectrometer1/Tb','Dimensions',{'channel_idx',retrievalTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTSys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/meanAngleAntenna','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/TRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

% Flags and flags
nccreate(filename,'/flags/time','Dimensions',{'time',Inf},'Datatype','double')
nccreate(filename,'/flags/flags','Dimensions',{'flags',length(calibratedSpectra(1).errorVector)},'Datatype','int64')

% We input a (1xerrorVectorSize) int vector to identify the errors
nccreate(filename,'/flags/calibration_flags','Dimensions',{'flags',length(calibratedSpectra(1).errorVector),'time',Inf},'Datatype','int64','FillValue',-9999)
%nccreate(filename,'/error/errorLabel','Dimensions',{'errorLabel',20,'time',Inf},'Datatype','char','FillValue','')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the spectrometer1 variable
ncwrite(filename,'/spectrometer1/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/spectrometer1/time','units',calibratedSpectra(1).meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time','calendar',calibratedSpectra(1).calendar);
ncwriteatt(filename,'/spectrometer1/time','description','mean time of the measurements for this cycle');

ncwrite(filename,'/spectrometer1/channel_idx',1:retrievalTool.numberOfChannels);

ncwrite(filename,'/spectrometer1/year',int64([calibratedSpectra.year]));
ncwrite(filename,'/spectrometer1/month',int64([calibratedSpectra.month]));
ncwrite(filename,'/spectrometer1/day',int64([calibratedSpectra.day]));
ncwrite(filename,'/spectrometer1/timeOfDay',[calibratedSpectra.timeOfDay]);

ncwrite(filename,'/spectrometer1/lat',retrievalTool.lat);
ncwrite(filename,'/spectrometer1/lon',retrievalTool.lon);
ncwrite(filename,'/spectrometer1/alt',retrievalTool.altitude);
%ncwrite(filename,'/spectrometer1/azimuthAngle', ??);

%ncwrite(filename,'/spectrometer1/startTime',char(startTimeInt8)');
%ncwrite(filename,'/spectrometer1/stopTime',datestr(calibratedSpectra(t).dateStop,'yyyymmddTHHMMSSZ'));

ncwrite(filename,'/spectrometer1/Tb',Tb');
ncwrite(filename,'/spectrometer1/THot',[calibratedSpectra.THot]);
ncwrite(filename,'/spectrometer1/stdTHot',[calibratedSpectra.stdTHot]);
ncwrite(filename,'/spectrometer1/TSys',[calibratedSpectra.Tsys]);
ncwrite(filename,'/spectrometer1/stdTSys',[calibratedSpectra.stdTSys]);
ncwrite(filename,'/spectrometer1/calibrationTime',[calibratedSpectra.calibrationTime]);
ncwrite(filename,'/spectrometer1/meanAngleAntenna',[calibratedSpectra.meanAngleAntenna]);

ncwrite(filename,'/spectrometer1/TRoom',[calibratedSpectra.TempRoom]);
ncwrite(filename,'/spectrometer1/TOut',[calibratedSpectra.TempOut]);

% Writing the errors variables
ncwrite(filename,'/flags/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/flags/time','units',calibratedSpectra(1).meanDatetimeUnit);
ncwriteatt(filename,'/flags/time','calendar',calibratedSpectra(1).calendar);
ncwriteatt(filename,'/flags/time','description','mean time of the measurements for this cycle');

ncwrite(filename,'/flags/flags',1:length(calibratedSpectra(1).errorVector));
ncwrite(filename,'/flags/calibration_flags',errorCalib');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the GLOBAL attributes of the files
% Originator attributes
ncwriteatt(filename,'/','title','Brightness temperature measured by ground-based radiometer');
ncwriteatt(filename,'/','location',retrievalTool.dataLocation);
ncwriteatt(filename,'/','source',retrievalTool.dataSource);
ncwriteatt(filename,'/','name',retrievalTool.PI_NAME);
ncwriteatt(filename,'/','institution',retrievalTool.PI_AFFILIATION);
ncwriteatt(filename,'/','contact',retrievalTool.PI_ADDRESS);
ncwriteatt(filename,'/','mail',retrievalTool.PI_EMAIL);
% ncwriteatt(filename,'/','DO_NAME','');
% ncwriteatt(filename,'/','DO_AFFILIATION','');
% ncwriteatt(filename,'/','DO_ADDRESS','');
% ncwriteatt(filename,'/','DO_EMAIL','');
% ncwriteatt(filename,'/','DS_NAME','');
% ncwriteatt(filename,'/','DS_AFFILIATION','');
% ncwriteatt(filename,'/','DS_ADDRESS','');
% ncwriteatt(filename,'/','DS_EMAIL','');
% Dataset attributes
%ncwriteatt(filename,'/','DATA_DISCIPLINE','ATMOSPHERIC.PHYSICS');
%ncwriteatt(filename,'/','DATA_GROUP','EXPERIMENTAL;...');

ncwriteatt(filename,'/','history','');
ncwriteatt(filename,'/','references','');
ncwriteatt(filename,'/','comment','');
%ncwriteatt(filename,'/','DATA_VARIABLES','');

ncwriteatt(filename,'/','rawData',log.file);
ncwriteatt(filename,'/','raw_data_software_version',num2str(log.SW_version(1)));
ncwriteatt(filename,'/','calibrated_version',calibratedSpectra(1).calibrationVersion);
ncwriteatt(filename,'/','raw_file_comment',log.comment);

% Geolocation attributes
ncwriteatt(filename,'/','data_start_date',datestr(calibratedSpectra(1).dateStart,'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','data_stop_date',datestr(calibratedSpectra(end).dateStop,'yyyymmddTHHMMSSZ'));
%ncwriteatt(filename,'/','DATA_FILE_VERSION','');
%ncwriteatt(filename,'/','DATA_RULES_OF_USE','');

% Create global variables --> saved as variables
%ncwriteatt(filename,'/','LATITUDE.INSTRUMENT',retrievalTool.lat);
%ncwriteatt(filename,'/','LONGITUDE.INSTRUMENT',retrievalTool.lon);
%ncwriteatt(filename,'/','ALTITUDE.INSTRUMENT',retrievalTool.altitude);

% Global file attributes
ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(now,'yyyymmddTHHMMSSZ'));   %TODO

ncwriteatt(filename,'/','rawFilename','')

ncwriteatt(filename,'/','featureType','timeSeries');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the attributes of the groups
% flags group
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables attributes
% The following are required for each variable
% Attribute name (GEOMS convention)
%attrName={'VAR_NAME','VAR_DESCRIPTION','VAR_SIZE','VAR_DEPEND','VAR_DATA_TYPE','VAR_UNITS','VAR_SI_CONVERSION','VAR_VALID_MIN','VAR_VALID_MAX','VAR_FILL_VALUE'};

% Attribute name (CF convention)
attrName={'long_name','standard_name','units','description'};

%Attributes for the spectrometer1 variables (CF convention)
attrVal.tod = {'TOD',...
    '',...
    'day',...
    'Time of the day'};

attrVal.lat = {'station latitude',...
    'latitude',...
    'degree_north',...
    'latitude defined according to WGS84'};

attrVal.lon = {'station longitude',...
    'longitude',...
    'degree_north',...
    'longitude defined according to WGS84'};

attrVal.alt = {'station altitude',...
    'altitude',...
    'm',...
    'above see level'};

attrVal.azimuth = {'azimuth angle',...
    'sensor_zenith_angle',...
    'degree',...
    'reference direction is NORTH'};

attrVal.Tb = {'Tb',...
    'brightness_temperature',...
    'K',...
    'calibrated brightness temperature for this cycle'};

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
    'min',...
    'Time used for calibrating the spectra'};

attrVal.meanAngleAntenna = {'meanAngleAntenna',...
    'elevation_angle',...
    'degree',...
    'mean elevation angle of the antenna during this cycle'};

attrVal.TRoom = {'TRoom',...
    '',...
    'K',...
    'mean room temperature'};

attrVal.TOut = {'TOut',...
    '',...
    'K',...
    'mean outside temperaure'};

for i=1:length(attrName)
    ncwriteatt(filename,'/spectrometer1/timeOfDay',attrName{i},attrVal.tod{i});
    ncwriteatt(filename,'/spectrometer1/lat',attrName{i},attrVal.lat{i});
    ncwriteatt(filename,'/spectrometer1/lon',attrName{i},attrVal.lon{i});
    ncwriteatt(filename,'/spectrometer1/alt',attrName{i},attrVal.alt{i});
    ncwriteatt(filename,'/spectrometer1/azimuthAngle',attrName{i},attrVal.azimuth{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/TSys',attrName{i},attrVal.TSys{i});
    ncwriteatt(filename,'/spectrometer1/stdTSys',attrName{i},attrVal.stdTSys{i});
    ncwriteatt(filename,'/spectrometer1/calibrationTime',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/meanAngleAntenna',attrName{i},attrVal.meanAngleAntenna{i});
    ncwriteatt(filename,'/spectrometer1/TRoom',attrName{i},attrVal.TRoom{i});
    ncwriteatt(filename,'/spectrometer1/TOut',attrName{i},attrVal.TOut{i});
end

if retrievalTool.saveAllCycles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding debug (optionnal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize matrices
    nClean=0;
    for t = 1:length(calibratedSpectra)
        nClean=nClean+calibratedSpectra(t).numberOfCleanAntennaAngle;
    end
    
    channel_idx=int64(ones(nClean,retrievalTool.numberOfChannels)*NaN);
    TbAll=single(ones(nClean,retrievalTool.numberOfChannels)*NaN);
    cycleNumber=int64(ones(nClean,1)*NaN);
    cycleId=int64(1:nClean);
    
    counter=1;
    for t = 1:length(calibratedSpectra)
        for i =1:calibratedSpectra(t).numberOfCleanAntennaAngle
            TbAll(counter,:)=calibratedSpectra(t).TbAll(i,:);
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

feedback=0;

end
