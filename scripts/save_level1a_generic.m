function feedback = save_level1a_generic(retrievalTool,log,calibratedSpectra,dateStr)
%
%
%
% Saving level1a into netCDF file
% We are using the netcdf package because it offers far more flexibility
% for the writing.
locationLevel1a=retrievalTool.level1Folder;

% Some additional variables to save in every netCDF files:
softwareVersion=log.SW_version(1);
rawDataName=log.file;
commentRawFile=log.comment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOK AT GEOMS CONVENTIONS 

% Here LOOOP in all calibration cycle
for t = 1:length(calibratedSpectra)
    filename=[locationLevel1a retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(t) '.h5'];
    %filename=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(t) '.h5'];
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
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create global variables
    h5create(filename,'/LATITUDE.INSTRUMENT',size(retrievalTool.altitude),'Datatype','double');
    h5create(filename,'/LONGITUDE.INSTRUMENT',size(retrievalTool.lon),'Datatype','double');
    h5create(filename,'/ALTITUDE.INSTRUMENT',size(retrievalTool.lat),'Datatype','double');
    h5create(filename,'/DATETIME',size(calibratedSpectra(t).meanDatetime),'Datatype','double');
    h5create(filename,'/YEAR',size(calibratedSpectra(t).year),'Datatype','double');
    h5create(filename,'/MONTH',size(calibratedSpectra(t).month),'Datatype','double');
    h5create(filename,'/DAY',size(calibratedSpectra(t).day),'Datatype','double');
    h5create(filename,'/TIMEOFDAY',size(calibratedSpectra(t).timeOfDay),'Datatype','double');
    
    % Create the different dataset
    % Scientific Dataset (SDS)
    h5create(filename,'/SDS/BRIGTHNESS_TEMPERATURE',size(calibratedSpectra(t).Tb),'Datatype','double')
    h5create(filename,'/SDS/MEAN_THOT',size(calibratedSpectra(t).THot),'Datatype','double')
    h5create(filename,'/SDS/STD_THOT',size(calibratedSpectra(t).stdTHot),'Datatype','double')
    h5create(filename,'/SDS/MEAN_TSYS',size(calibratedSpectra(t).Tsys),'Datatype','double')
    h5create(filename,'/SDS/STD_TSYS',size(calibratedSpectra(t).stdTSys),'Datatype','double')
    h5create(filename,'/SDS/CALIBRATION_TIME',size(calibratedSpectra(t).calibrationTime),'Datatype','double')
    h5create(filename,'/SDS/MEAN_ANGLE_ANTENNA',size(calibratedSpectra(t).meanAngleAntenna),'Datatype','double')
    
    h5create(filename,'/SDS/MEAN_TROOM',size(calibratedSpectra(t).TempRoom),'Datatype','double')
    h5create(filename,'/SDS/MEAN_TOUT',size(calibratedSpectra(t).TempOut),'Datatype','double')
    
    % Flags and error
    h5create(filename,'/ERROR/calibrationError',1,'Datatype','double')
    
    % Writing the GLOBAL attributes of the files
    %varid = netcdf.getConstant('NC_GLOBAL');
    % Originator attributes
    h5writeatt(filename,'/','PI_NAME',retrievalTool.PI_NAME);
    h5writeatt(filename,'/','PI_AFFILIATION',retrievalTool.PI_AFFILIATION);
    h5writeatt(filename,'/','PI_ADDRESS',retrievalTool.PI_ADDRESS);
    h5writeatt(filename,'/','PI_EMAIL',retrievalTool.PI_EMAIL);
    h5writeatt(filename,'/','DO_NAME','');
    h5writeatt(filename,'/','DO_AFFILIATION','');
    h5writeatt(filename,'/','DO_ADDRESS','');
    h5writeatt(filename,'/','DO_EMAIL','');
    h5writeatt(filename,'/','DS_NAME','');
    h5writeatt(filename,'/','DS_AFFILIATION','');
    h5writeatt(filename,'/','DS_ADDRESS','');
    h5writeatt(filename,'/','DS_EMAIL','');
    % Dataset attributes
    h5writeatt(filename,'/','DATA_DESCRIPTION','Brightness temperature measured by ground-based radiometer');
    h5writeatt(filename,'/','DATA_DISCIPLINE','ATMOSPHERIC.PHYSICS');
    h5writeatt(filename,'/','DATA_GROUP','EXPERIMENTAL;...');
    h5writeatt(filename,'/','DATA_LOCATION',retrievalTool.dataLocation);
    h5writeatt(filename,'/','DATA_SOURCE',retrievalTool.dataSource);
    h5writeatt(filename,'/','DATA_VARIABLES','');
    
    h5writeatt(filename,'/','RAW_DATA',log.file);
    h5writeatt(filename,'/','RAW_DATA_SOFTWARE_VERSION',num2str(log.SW_version(1)));
    h5writeatt(filename,'/','CALIBRATION_VERSION',calibratedSpectra(1).calibrationVersion);
    h5writeatt(filename,'/','RAW_FILE_COMMENT',log.comment);
    
    % Geolocation attributes
    h5writeatt(filename,'/','DATA_START_DATE',datestr(calibratedSpectra(t).dateStart,'yyyymmddTHHMMSSZ'));
    h5writeatt(filename,'/','DATA_STOP_DATE',datestr(calibratedSpectra(t).dateStop,'yyyymmddTHHMMSSZ'));
    h5writeatt(filename,'/','DATA_FILE_VERSION','');
    h5writeatt(filename,'/','DATA_RULES_OF_USE','');
    
    % Global file attributes
    h5writeatt(filename,'/','FILE_NAME',filename);
    h5writeatt(filename,'/','FILE_GENERATION_DATE',datestr(now,'yyyymmddTHHMMSSZ'));   %TODO
    
    h5writeatt(filename,'/','RAW_FILENAME','');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writing hdf5 variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Writing the GLOBAL variables
    h5write(filename,'/LATITUDE.INSTRUMENT',retrievalTool.lat)
    h5write(filename,'/LONGITUDE.INSTRUMENT',retrievalTool.lon)
    h5write(filename,'/ALTITUDE.INSTRUMENT',retrievalTool.altitude);
    h5write(filename,'/DATETIME',calibratedSpectra(t).meanDatetime);
    h5write(filename,'/YEAR',calibratedSpectra(t).year);
    h5write(filename,'/MONTH',calibratedSpectra(t).month);
    h5write(filename,'/DAY',calibratedSpectra(t).day);
    h5write(filename,'/TIMEOFDAY',calibratedSpectra(t).timeOfDay);
    
    % Writing the SDS variable
    h5write(filename,'/SDS/BRIGTHNESS_TEMPERATURE',calibratedSpectra(t).Tb);
    h5write(filename,'/SDS/MEAN_THOT',calibratedSpectra(t).THot);
    h5write(filename,'/SDS/STD_THOT',calibratedSpectra(t).stdTHot);
    h5write(filename,'/SDS/MEAN_TSYS',calibratedSpectra(t).Tsys);
    h5write(filename,'/SDS/STD_TSYS',calibratedSpectra(t).stdTSys);
    h5write(filename,'/SDS/CALIBRATION_TIME',calibratedSpectra(t).calibrationTime);
    h5write(filename,'/SDS/MEAN_ANGLE_ANTENNA',calibratedSpectra(t).meanAngleAntenna);
    
    h5write(filename,'/SDS/MEAN_TROOM',calibratedSpectra(t).TempRoom)
    h5write(filename,'/SDS/MEAN_TOUT',calibratedSpectra(t).TempOut)
    
    % Writing the errors variables
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Variables attributes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following are required for each variable
    % VAR_NAME x cf. Section 5.1.1
    % VAR_DESCRIPTION x cf. Section 5.1.2
    % VAR_NOTES o cf. Section 5.1.3
    % VAR_SIZE x cf. Section 5.1.4
    % VAR_DEPEND x cf. Section 5.1.5
    % VAR_DATA_TYPE x http://avdc.gsfc.nasa.gov/GEOMS.php?GEOMS=VAR_DATA_TYPE
    % VAR_UNITS x http://avdc.gsfc.nasa.gov/GEOMS.php?GEOMS=VAR_UNITS
    % VAR_SI_CONVERSION x cf. Section 5.1.8
    % VAR_VALID_MIN x cf. Section 5.1.9
    % VAR_VALID_MAX x cf. Section 5.1.10
    % VAR_FILL_VALUE
    
    % Attribute name:
    attrName={'VAR_NAME','VAR_DESCRIPTION','VAR_SIZE','VAR_DEPEND','VAR_DATA_TYPE','VAR_UNITS','VAR_SI_CONVERSION','VAR_VALID_MIN','VAR_VALID_MAX','VAR_FILL_VALUE'};
    
    % Attributes for the GLOBAL variables   
    attrVal.Latitude = {'LATITUDE.INSTRUMENT',...
        'Latitude of the observation site (deg)',...
        '1',...
        'CONSTANT',...
        'double',...
        'deg',...
        '0.0;1.74533E-2;rad', ...
        '-90',...
        '90',...
        '-9999'};
    
    attrVal.Longitude = {'LONGITUDE.INSTRUMENT',...
        'Longitude of the observation site (deg)',...
        '1',...
        'CONSTANT',...
        'double',...
        'deg',...
        '0.0;1.74533E-2;rad', ...
        '-180',...
        '180',...
        '-9999'};
    
    % Attributes for the SDS variables
    % Tb
    attrVal.Tb = {'Tb',...
        'Brightness Temperature',...
        '1',...
        'DATETIME',...
        '',...
        'K',...
        '0.0;1.0;K', ...
        '0',...
        '1000000',...
        '-9999'};
    
    for i=1:length(attrName)
        h5writeatt(filename,'/LATITUDE.INSTRUMENT',attrName{i},attrVal.Latitude{i});
        h5writeatt(filename,'/LONGITUDE.INSTRUMENT',attrName{i},attrVal.Longitude{i});
        h5writeatt(filename,'/SDS/BRIGTHNESS_TEMPERATURE',attrName{i},attrVal.Tb{i});
    end
    
    %h5writeatt(filename,'/SDS/BRIGTHNESS_TEMPERATURE','VAR_UNITS','K');
    
    
end
%h5writeatt(filename,'/FFT/Tb','unit','Kelvin');

% 
% h5writeatt(filename,'/','source','Ground-based radiometer');
% h5writeatt(filename,'/','instrument',retrievalTool.instrumentName);
% 
% h5writeatt(filename,'/','history',);   %TODO
% h5writeatt(filename,'/','references','');
% h5writeatt(filename,'/','comment','WORK IN PROGRESS');
% h5writeatt(filename,'/','conventions','...');

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
