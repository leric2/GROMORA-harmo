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

%initialize some variable:
meanTime=ones(length(calibratedSpectra),1)*NaN;
channelId=ones(length(calibratedSpectra),retrievalTool.numberOfChannels)*NaN;
Tb=ones(length(calibratedSpectra),retrievalTool.numberOfChannels)*NaN;

error=ones(length(calibratedSpectra),length(calibratedSpectra(1).errorVector))*NaN;
errorCalib=ones(length(calibratedSpectra),length(calibratedSpectra(1).errorVector))*NaN;

for t = 1:length(calibratedSpectra)
    meanTime(t)=calibratedSpectra(t).meanDatetime;
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
nccreate(filename,'/spectrometer1/channel_idx','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/year','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/month','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/day','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/timeOfDay','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
%nccreate(filename,'/spectrometer1/startTime','Dimensions',{'time',Inf},'Datatype','char','FillValue','NA');
%nccreate(filename,'/spectrometer1/stopTime','Dimensions',{'time',Inf},'Datatype','char','FillValue','NA');

nccreate(filename,'/spectrometer1/Tb','Dimensions',{'channel_idx',retrievalTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/meanTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/meanTsys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTsys','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/meanAngleAntenna','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/meanTRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/meanTOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

% Flags and error
nccreate(filename,'/error/time','Dimensions',{'time',Inf},'Datatype','double')
nccreate(filename,'/error/error','Dimensions',{'error',length(calibratedSpectra(1).errorVector)},'Datatype','double')

% We input a (1xerrorVectorSize) int vector to identify the errors
nccreate(filename,'/error/calibration_error','Dimensions',{'error',length(calibratedSpectra(1).errorVector),'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the spectrometer1 variable
ncwrite(filename,'/spectrometer1/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/spectrometer1/time','units',calibratedSpectra(1).meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time','calendar',calibratedSpectra(1).calendar);

ncwrite(filename,'/spectrometer1/channel_idx',1:retrievalTool.numberOfChannels);

ncwrite(filename,'/spectrometer1/year',[calibratedSpectra.year]);
ncwrite(filename,'/spectrometer1/month',[calibratedSpectra.month]);
ncwrite(filename,'/spectrometer1/day',[calibratedSpectra.day]);
ncwrite(filename,'/spectrometer1/timeOfDay',[calibratedSpectra.timeOfDay]);

%ncwrite(filename,'/spectrometer1/startTime',datestr(calibratedSpectra(t).dateStart,'yyyymmddTHHMMSSZ'));
%ncwrite(filename,'/spectrometer1/stopTime',datestr(calibratedSpectra(t).dateStop,'yyyymmddTHHMMSSZ'));

ncwrite(filename,'/spectrometer1/Tb',Tb');
ncwrite(filename,'/spectrometer1/meanTHot',[calibratedSpectra.THot]);
ncwrite(filename,'/spectrometer1/stdTHot',[calibratedSpectra.stdTHot]);
ncwrite(filename,'/spectrometer1/meanTsys',[calibratedSpectra.Tsys]);
ncwrite(filename,'/spectrometer1/stdTsys',[calibratedSpectra.stdTSys]);
ncwrite(filename,'/spectrometer1/calibrationTime',[calibratedSpectra.calibrationTime]);
ncwrite(filename,'/spectrometer1/meanAngleAntenna',[calibratedSpectra.meanAngleAntenna]);

ncwrite(filename,'/spectrometer1/meanTRoom',[calibratedSpectra.TempRoom]);
ncwrite(filename,'/spectrometer1/meanTOut',[calibratedSpectra.TempOut]);

% Writing the errors variables
ncwrite(filename,'/error/time',meanTime);
ncwriteatt(filename,'/error/time','units',calibratedSpectra(1).meanDatetimeUnit);
ncwriteatt(filename,'/error/time','calendar',calibratedSpectra(1).calendar);

ncwrite(filename,'/error/error',1:length(calibratedSpectra(1).errorVector));
ncwrite(filename,'/error/calibration_error',errorCalib');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the GLOBAL attributes of the files
%varid = netcdf.getConstant('NC_GLOBAL');
% Originator attributes
ncwriteatt(filename,'/','PI_NAME',retrievalTool.PI_NAME);
ncwriteatt(filename,'/','PI_AFFILIATION',retrievalTool.PI_AFFILIATION);
ncwriteatt(filename,'/','PI_ADDRESS',retrievalTool.PI_ADDRESS);
ncwriteatt(filename,'/','PI_EMAIL',retrievalTool.PI_EMAIL);
ncwriteatt(filename,'/','DO_NAME','');
ncwriteatt(filename,'/','DO_AFFILIATION','');
ncwriteatt(filename,'/','DO_ADDRESS','');
ncwriteatt(filename,'/','DO_EMAIL','');
ncwriteatt(filename,'/','DS_NAME','');
ncwriteatt(filename,'/','DS_AFFILIATION','');
ncwriteatt(filename,'/','DS_ADDRESS','');
ncwriteatt(filename,'/','DS_EMAIL','');
% Dataset attributes
ncwriteatt(filename,'/','DATA_DESCRIPTION','Brightness temperature measured by ground-based radiometer');
ncwriteatt(filename,'/','DATA_DISCIPLINE','ATMOSPHERIC.PHYSICS');
ncwriteatt(filename,'/','DATA_GROUP','EXPERIMENTAL;...');
ncwriteatt(filename,'/','DATA_LOCATION',retrievalTool.dataLocation);
ncwriteatt(filename,'/','DATA_SOURCE',retrievalTool.dataSource);
ncwriteatt(filename,'/','DATA_VARIABLES','');

ncwriteatt(filename,'/','RAW_DATA',log.file);
ncwriteatt(filename,'/','RAW_DATA_SOFTWARE_VERSION',num2str(log.SW_version(1)));
ncwriteatt(filename,'/','CALIBRATION_VERSION',calibratedSpectra(1).calibrationVersion);
ncwriteatt(filename,'/','RAW_FILE_COMMENT',log.comment);

% Geolocation attributes
ncwriteatt(filename,'/','DATA_START_DATE',datestr(calibratedSpectra(1).dateStart,'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','DATA_STOP_DATE',datestr(calibratedSpectra(end).dateStop,'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','DATA_FILE_VERSION','');
ncwriteatt(filename,'/','DATA_RULES_OF_USE','');

% Create global variables
ncwriteatt(filename,'/','LATITUDE.INSTRUMENT',retrievalTool.lat);
ncwriteatt(filename,'/','LONGITUDE.INSTRUMENT',retrievalTool.lon);
ncwriteatt(filename,'/','ALTITUDE.INSTRUMENT',retrievalTool.altitude);

% Global file attributes
ncwriteatt(filename,'/','FILE_NAME',filename);
ncwriteatt(filename,'/','FILE_GENERATION_DATE',datestr(now,'yyyymmddTHHMMSSZ'));   %TODO

ncwriteatt(filename,'/','RAW_FILENAME','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the attributes of the groups
% error group
ncwriteatt(filename,'/error','sufficientNumberOfIndices','');
ncwriteatt(filename,'/error','systemTemperatureOK','');
ncwriteatt(filename,'/error','hotAngleRemoved','');
ncwriteatt(filename,'/error','coldAngleRemoved','');
ncwriteatt(filename,'/error','antennaAngleRemoved','');
ncwriteatt(filename,'/error','LN2SensorsOK','');
ncwriteatt(filename,'/error','LN2LevelOK','');
ncwriteatt(filename,'/error','hotLoadOK','');
ncwriteatt(filename,'/error','FFT_adc_overload_OK','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables attributes
% The following are required for each variable
% Attribute name (GEOMS convention)
%attrName={'VAR_NAME','VAR_DESCRIPTION','VAR_SIZE','VAR_DEPEND','VAR_DATA_TYPE','VAR_UNITS','VAR_SI_CONVERSION','VAR_VALID_MIN','VAR_VALID_MAX','VAR_FILL_VALUE'};

% Attribute name (CF convention)
attrName={'standard_name','long_name','units','description'};

%Attributes for the spectrometer1 variables (CF convention)
% Tb
attrVal.Tb = {'Tb',...
    'Brightness Temperature',...
    'K',...
    ''};

attrVal.meanThot = {'meanTHot',...
    'Temperature of the hot load',...
    'K',...
    ''};

for i=1:length(attrName)
    %ncwriteatt(filename,'/LATITUDE.INSTRUMENT',attrName{i},attrVal.Latitude{i});
    ncwriteatt(filename,'/spectrometer1/Tb',attrName{i},attrVal.Tb{i});
    ncwriteatt(filename,'/spectrometer1/meanTHot',attrName{i},attrVal.Tb{i});
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
    
    channel_idx=ones(nClean,retrievalTool.numberOfChannels)*NaN;
    TbAll=ones(nClean,retrievalTool.numberOfChannels)*NaN;
    cycleNumber=ones(nClean,1)*NaN;
    cycleId=1:nClean;
    
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
    nccreate(filename,'/debug/cycle_id','Dimensions',{'cycle_id',Inf},'Datatype','double')
    nccreate(filename,'/debug/channel_idx','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    
    nccreate(filename,'/debug/calibration_cycle_number','Dimensions',{'cycle_id',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/debug/Tb_all','Dimensions',{'channel_idx',retrievalTool.numberOfChannels,'cycle_id',Inf},'Datatype','double','FillValue',-9999)
    
    
    % Writing the debug variables
    ncwrite(filename,'/debug/cycle_id',cycleId);
    
    ncwrite(filename,'/debug/channel_idx',1:retrievalTool.numberOfChannels);
    ncwrite(filename,'/debug/calibration_cycle_number',cycleNumber');
    ncwrite(filename,'/debug/Tb_all',TbAll');
end

feedback=0;

end
