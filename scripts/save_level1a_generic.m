function save_level1a_generic(retrievalTool,log,calibratedSpectra,dateStr)
%
%
%
% Saving level1a into netCDF file
% We are using the netcdf package because it offers far more flexibility
% for the writing.

calibrationVersion=0; % TODO as input of this function

locationLevel1a=retrievalTool.level1Folder;


% Some additional variables to save in every netCDF files:
softwareVersion=log.SW_version(1);
rawDataName=log.file;
commentRawFile=log.comment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOK AT GEOMS CONVENTIONS 

% Here LOOOP in all calibration cycle
% for i = 1:nCalir...............
i=1;

yyyy = calibratedSpectra(i).dateStart(1);
mm = calibratedSpectra(i).dateStart(2);
dd = calibratedSpectra(i).dateStart(3);
HH = calibratedSpectra(i).dateStart(4);
MM = calibratedSpectra(i).dateStart(5);
SS = calibratedSpectra(i).dateStart(6);

filename=[locationLevel1a retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(i) '.h5'];
title=[retrievalTool.instrumentName '_level1a_' dateStr '_' num2str(i)];

if exist(filename)
    %ncid=netcdf.open(filename,'WRITE');
    delete(filename)
else
    % Define the creation mode of our files
    %cmode = netcdf.getConstant('NETCDF4');
    %cmode = bitor(cmode,netcdf.getConstant('CLOBBER'));
    
    % Create the files
    %ncid = netcdf.create(filename,cmode);
end
% nccreate(filename,'Tb','Dimensions',{'Brightness Temperature',1,'Channels',retrievalTool.numberOfChannels},'Format','netcdf4')

% ncwrite(filename,'Tb',calibratedSpectra(1).Tb)

% Create some global variable
h5create(filename,'/lat',1,'Datatype','double')
h5create(filename,'/lon',1,'Datatype','double')
h5create(filename,'/time',1,'Datatype','double')

% Writing the GLOBAL attributes of the files
%varid = netcdf.getConstant('NC_GLOBAL');
h5writeatt(filename,'/','title',title);
h5writeatt(filename,'/','institution','Institude of Applied Physics, Universtiy of Bern, Switzerland');
h5writeatt(filename,'/','contact','eric.sauveat@iap.unibe.ch');
h5writeatt(filename,'/','source','Ground-based radiometer');
h5writeatt(filename,'/','instrument',retrievalTool.instrumentName);
h5writeatt(filename,'/','site',retrievalTool.siteName);
h5writeatt(filename,'/','history',datestr(now));   %TODO
h5writeatt(filename,'/','references','');
h5writeatt(filename,'/','comment','WORK IN PROGRESS');
h5writeatt(filename,'/','conventions','...');

% Create the different dataset
% Data
h5create(filename,'/FFT/Tb',size(calibratedSpectra(i).Tb),'Datatype','double')
h5create(filename,'/FFT/meanTsys',size(calibratedSpectra(i).Tb),'Datatype','double')

% Flag
h5create(filename,'/error/calibrationError',1,'Datatype','double')


% Writing some GLOBAL variable
% Latitude of the observation site
%dimid=netcdf.defDim(ncid,'lat',1);
%varID=netcdf.defVar(ncid,'lat','NC_FLOAT',dimid);
%netcdf.putVar(ncid,varID,retrievalTool.lat);
h5write(filename,'/lat',retrievalTool.lat)
h5write(filename,'/lon',retrievalTool.lon)

% Writing the group variable
h5write(filename,'/FFT/Tb',calibratedSpectra(i).Tb)


% Attributes for the variables
h5writeatt(filename,'/FFT/Tb','unit','Kelvin');

%netcdf.putAtt(ncid,varid,'site_lon',retrievalTool.lon);
%netcdf.putAtt(ncid,varid,'altitude',retrievalTool.altitude);

% Creating 2 groups in the file, one for the data, one for the flags
%spectroGroupID = netcdf.defGrp(ncid,'spectrometer0');
%errorGroupID = netcdf.defGrp(ncid,'level0-1a errors');


end
