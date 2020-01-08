function save_level1a_generic(retrievalTool,calibratedSpectra,dateStr)
%
%
%
% Saving level1a into netCDF file
% We are using the netcdf package because it offers far more flexibility
% for the writing.


locationLevel1a=retrievalTool.level1Folder;

filename=[locationLevel1a retrievalTool.instrumentName '_level1a_' dateStr '.nc'];

if exist(filename)
    ncid=netcdf.open(filename,'WRITE');
    
else
    % Define the creation mode of our files
    cmode = netcdf.getConstant('NETCDF4');
    cmode = bitor(cmode,netcdf.getConstant('NOCLOBBER'));
    
    % Create the files
    ncid = netcdf.create(filename,cmode);
end
% nccreate(filename,'Tb','Dimensions',{'Brightness Temperature',1,'Channels',retrievalTool.numberOfChannels},'Format','netcdf4')

% ncwrite(filename,'Tb',calibratedSpectra(1).Tb)

title=[retrievalTool.instrumentName '_level1a_' dateStr];

% Writing the GLOBAL attributes of the files
varid = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,varid,'title',title);
netcdf.putAtt(ncid,varid,'institution','Institude of Applied Physics, Universtiy of Bern, Switzerland');
netcdf.putAtt(ncid,varid,'contact','eric.sauveat@iap.unibe.ch');
netcdf.putAtt(ncid,varid,'source','Ground-based radiometer');
netcdf.putAtt(ncid,varid,'instrument',retrievalTool.instrumentName);
netcdf.putAtt(ncid,varid,'site',retrievalTool.siteName);

netcdf.putAtt(ncid,varid,'history',datestr(now));   %TODO
netcdf.putAtt(ncid,varid,'references','');
netcdf.putAtt(ncid,varid,'comment','WORK IN PROGRESS');
netcdf.putAtt(ncid,varid,'Conventions','CF-1.8');

% Writing some GLOBAL variable
% Latitude of the observation site
dimid=netcdf.defDim(ncid,'lat',1);
netcdf.defVar(ncid,'lat','NC_FLOAT',dimid);

netcdf.putAtt(ncid,varid,'site_lon',retrievalTool.lon);
netcdf.putAtt(ncid,varid,'altitude',retrievalTool.altitude);

% Creating 2 groups in the file, one for the data, one for the flags
spectroGroupID = netcdf.defGrp(ncid,'spectrometer0');
errorGroupID = netcdf.defGrp(ncid,'level0-1a errors');

netcdf.close(ncid);
end
