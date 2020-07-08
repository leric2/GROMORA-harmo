function feedback = save_level1a_generic(calibrationTool,log,calibratedSpectra)
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
% Saving level1a into netCDF file
% We are using the netcdf package because it offers far more flexibility
% for the writing.
locationLevel1a=calibrationTool.level1Folder;

% Some additional variables to save in every netCDF files:
softwareVersion=log.SW_version(1);
rawDataName=log.file;
commentRawFile=log.comment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOK AT GEOMS CONVENTIONS 

% Here LOOOP in all calibration cycle
for t = 1:length(calibratedSpectra)
    filename=[locationLevel1a calibrationTool.instrumentName '_level1a_' calibrationTool.dateStr '_' sprintf('%02d',t) '.nc'];
    %filename=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(t) '.nc'];
    %title=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(i)];
    
    if exist(filename)
        %ncid=netcdf.open(filename,'WRITE');
        delete(filename)
        % Define the creation mode of our files
        %cmode = netcdf.getConstant('NETCDF4');
        %cmode = bitor(cmode,netcdf.getConstant('CLOBBER'));
        
        % Create the files
        %ncid = netcdf.create(filename,cmode);
    end
    % nccreate(filename,'Tb','Dimensions',{'Brightness Temperature',1,'Channels',retrievalTool.numberOfChannels},'Format','netcdf4')
    
    % ncwrite(filename,'Tb',calibratedSpectra(1).Tb)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Create the different dataset
    % Scientific Dataset (SDS)
    % First create coordinates variable (enable 'netcdf4' format)
    nccreate(filename,'/SDS/time','Dimensions',{'time',Inf},'Datatype','double','Format','netcdf4');
    nccreate(filename,'/SDS/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    
    nccreate(filename,'/SDS/year','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/SDS/month','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/SDS/day','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/SDS/timeOfDay','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
    %nccreate(filename,'/SDS/startTime','Dimensions',{'time',Inf},'Datatype','char','FillValue','NA');
    %nccreate(filename,'/SDS/stopTime','Dimensions',{'time',Inf},'Datatype','char','FillValue','NA');
    
    nccreate(filename,'/SDS/Tb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/SDS/meanTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/SDS/STD_THOT','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/SDS/MEAN_TSYS','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/SDS/STD_TSYS','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/SDS/CALIBRATION_TIME','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/SDS/MEAN_ANGLE_ANTENNA','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    
    nccreate(filename,'/SDS/MEAN_TROOM','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    nccreate(filename,'/SDS/MEAN_TOUT','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
    
    % Flags and error
    nccreate(filename,'/error/time','Dimensions',{'time',Inf},'Datatype','double')
    nccreate(filename,'/error/error','Dimensions',{'error',length(calibratedSpectra(t).errorVector)},'Datatype','double')
    
    % We input a (1xerrorVectorSize) int vector to identify the errors
    nccreate(filename,'/error/calibration_error','Dimensions',{'error',length(calibratedSpectra(t).errorVector),'time',Inf},'Datatype','double','FillValue',-9999)
    
    % Group for debugging variables:
    %nccreate(filename,'/debug/cycle_id','Dimensions',{'cycle_id',Inf},'Datatype','double')
    %nccreate(filename,'/debug/channel_idx','Dimensions',{'channel_idx',retrievalTool.numberOfChannels},'Datatype','double','FillValue',-9999)
    
    %nccreate(filename,'/debug/cycle_time','Dimensions',{'cycle_id',Inf},'Datatype','double','FillValue',-9999)
    
    % Writing the GLOBAL attributes of the files
    %varid = netcdf.getConstant('NC_GLOBAL');
    % Originator attributes
    ncwriteatt(filename,'/','PI_NAME',calibrationTool.PI_NAME);
    ncwriteatt(filename,'/','PI_AFFILIATION',calibrationTool.PI_AFFILIATION);
    ncwriteatt(filename,'/','PI_ADDRESS',calibrationTool.PI_ADDRESS);
    ncwriteatt(filename,'/','PI_EMAIL',calibrationTool.PI_EMAIL);
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
    ncwriteatt(filename,'/','DATA_LOCATION',calibrationTool.dataLocation);
    ncwriteatt(filename,'/','DATA_SOURCE',calibrationTool.dataSource);
    ncwriteatt(filename,'/','DATA_VARIABLES','');
    
    ncwriteatt(filename,'/','RAW_DATA',log.file);
    ncwriteatt(filename,'/','RAW_DATA_SOFTWARE_VERSION',num2str(log.SW_version(1)));
    ncwriteatt(filename,'/','CALIBRATION_VERSION',calibratedSpectra(1).calibrationVersion);
    ncwriteatt(filename,'/','RAW_FILE_COMMENT',log.comment);
    
    % Geolocation attributes
    ncwriteatt(filename,'/','DATA_START_DATE',datestr(calibratedSpectra(t).dateStart,'yyyymmddTHHMMSSZ'));
    ncwriteatt(filename,'/','DATA_STOP_DATE',datestr(calibratedSpectra(t).dateStop,'yyyymmddTHHMMSSZ'));
    ncwriteatt(filename,'/','DATA_FILE_VERSION','');
    ncwriteatt(filename,'/','DATA_RULES_OF_USE','');
    
    % Create global variables
    ncwriteatt(filename,'/','LATITUDE.INSTRUMENT',calibrationTool.lat);
    ncwriteatt(filename,'/','LONGITUDE.INSTRUMENT',calibrationTool.lon);
    ncwriteatt(filename,'/','ALTITUDE.INSTRUMENT',calibrationTool.altitude);
    
    % Global file attributes
    ncwriteatt(filename,'/','FILE_NAME',filename);
    ncwriteatt(filename,'/','FILE_GENERATION_DATE',datestr(now,'yyyymmddTHHMMSSZ'));   %TODO
    
    ncwriteatt(filename,'/','RAW_FILENAME','');
    
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writing netCDF variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writing the SDS variable
    ncwrite(filename,'/SDS/time',calibratedSpectra(t).meanDatetime);
    ncwriteatt(filename,'/SDS/time','units',calibratedSpectra(t).meanDatetimeUnit);
    ncwriteatt(filename,'/SDS/time','calendar',calibratedSpectra(t).calendar);
    
    ncwrite(filename,'/SDS/channel_idx',1:calibrationTool.numberOfChannels);
    
    ncwrite(filename,'/SDS/year',calibratedSpectra(t).year);
    ncwrite(filename,'/SDS/month',calibratedSpectra(t).month);
    ncwrite(filename,'/SDS/day',calibratedSpectra(t).day);
    ncwrite(filename,'/SDS/timeOfDay',calibratedSpectra(t).timeOfDay);
    
    %ncwrite(filename,'/SDS/startTime',datestr(calibratedSpectra(t).dateStart,'yyyymmddTHHMMSSZ'));
    %ncwrite(filename,'/SDS/stopTime',datestr(calibratedSpectra(t).dateStop,'yyyymmddTHHMMSSZ'));
    
    ncwrite(filename,'/SDS/Tb',calibratedSpectra(t).Tb');
    ncwrite(filename,'/SDS/meanTHot',calibratedSpectra(t).THot);
    ncwrite(filename,'/SDS/STD_THOT',calibratedSpectra(t).stdTHot);
    ncwrite(filename,'/SDS/MEAN_TSYS',calibratedSpectra(t).Tsys);
    ncwrite(filename,'/SDS/STD_TSYS',calibratedSpectra(t).stdTSys);
    ncwrite(filename,'/SDS/CALIBRATION_TIME',calibratedSpectra(t).calibrationTime);
    ncwrite(filename,'/SDS/MEAN_ANGLE_ANTENNA',calibratedSpectra(t).meanAngleAntenna);
    
    ncwrite(filename,'/SDS/MEAN_TROOM',calibratedSpectra(t).TempRoom)
    ncwrite(filename,'/SDS/MEAN_TOUT',calibratedSpectra(t).TempOut)
    
    % Writing the errors variables
    ncwrite(filename,'/error/time',calibratedSpectra(t).meanDatetime);
    ncwriteatt(filename,'/error/time','units',calibratedSpectra(t).meanDatetimeUnit);
    ncwriteatt(filename,'/error/time','calendar',calibratedSpectra(t).calendar);
    
    ncwrite(filename,'/error/error',(1:length(calibratedSpectra(t).errorVector)));
    ncwrite(filename,'/error/calibration_error',calibratedSpectra(t).errorVector');
    
     
    % Writing the debug variables
    %ncwrite(filename,'/debug/cycle_id',);
   
    %ncwrite(filename,'/debug/cycle_time',);
    %ncwrite(filename,'/debug/Tb_all',);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Variables attributes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following are required for each variable
    % Attribute name (GEOMS convention)
    %attrName={'VAR_NAME','VAR_DESCRIPTION','VAR_SIZE','VAR_DEPEND','VAR_DATA_TYPE','VAR_UNITS','VAR_SI_CONVERSION','VAR_VALID_MIN','VAR_VALID_MAX','VAR_FILL_VALUE'};
    
    % Attribute name (CF convention)
    attrName={'standard_name','long_name','units','description'};
    
    %Attributes for the SDS variables (CF convention)
    % Tb
    attrVal.Tb = {'Tb',...
        'Brightness Temperature',...
        'K',...
        ''};
    
    attrVal.meanThot = {'meanTHot',...
        'Temperature of the hot load',...
        'K',...
        ''};
    
%        
%     attrVal.Latitude = {'LATITUDE.INSTRUMENT',...
%         'Latitude of the observation site (deg)',...
%         '1',...
%         'CONSTANT',...
%         'double',...
%         'deg',...
%         '0.0;1.74533E-2;rad', ...
%         '-90',...
%         '90',...
%         '-9999'};
%     
%     attrVal.Longitude = {'LONGITUDE.INSTRUMENT',...
%         'Longitude of the observation site (deg)',...
%         '1',...
%         'CONSTANT',...
%         'double',...
%         'deg',...
%         '0.0;1.74533E-2;rad', ...
%         '-180',...
%         '180',...
%         '-9999'};
%     
%     % Attributes for the SDS variables
%     % Tb
%     attrVal.Tb = {'Tb',...
%         'Brightness Temperature',...
%         '1',...
%         'DATETIME',...
%         '',...
%         'K',...
%         '0.0;1.0;K', ...
%         '0',...
%         '1000000',...
%         '-9999'};
    %ncwriteatt(filename,'/LATITUDE.INSTRUMENT','_FillValue','NaN');
    
    for i=1:length(attrName)
        %ncwriteatt(filename,'/LATITUDE.INSTRUMENT',attrName{i},attrVal.Latitude{i});
        ncwriteatt(filename,'/SDS/Tb',attrName{i},attrVal.Tb{i});
        ncwriteatt(filename,'/SDS/meanTHot',attrName{i},attrVal.Tb{i});
    end
    
    %ncwriteatt(filename,'/SDS/BRIGTHNESS_TEMPERATURE','VAR_UNITS','K');
    
    
end
%ncwriteatt(filename,'/FFT/Tb','unit','Kelvin');

% 
% ncwriteatt(filename,'/','source','Ground-based radiometer');
% ncwriteatt(filename,'/','instrument',retrievalTool.instrumentName);
% 
% ncwriteatt(filename,'/','history',);   %TODO
% ncwriteatt(filename,'/','references','');
% ncwriteatt(filename,'/','comment','WORK IN PROGRESS');
% ncwriteatt(filename,'/','conventions','...');

%dimid=netcdf.defDim(ncid,'lat',1);
%varID=netcdf.defVar(ncid,'lat','NC_FLOAT',dimid);
%netcdf.putVar(ncid,varID,retrievalTool.lat);


%netcdf.putAtt(ncid,varid,'site_lon',retrievalTool.lon);
%netcdf.putAtt(ncid,varid,'altitude',retrievalTool.altitude);

% Creating 2 groups in the file, one for the data, one for the flags
%spectroGroupID = netcdf.defGrp(ncid,'spectrometer0');
%errorGroupID = netcdf.defGrp(ncid,'level0-1a errors');
feedback=0;

end
