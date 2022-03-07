function calibrationTool = save_level1a_daily(calibrationTool,logFile,calibratedSpectra,warningLevel0)
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
if calibrationTool.calType=="debug"
    filename=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '_debug.nc'];
else
    filename=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr calibrationTool.extraName '.nc'];
end
calibrationTool.filenameLevel1a=filename;

% Rewrite if already existing
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
% Calibration variables   
%nccreate(filename,'/spectrometer1/effectiveCalibrationTime','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/Tb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/stdTb','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
nccreate(filename,'/spectrometer1/frequencies','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/intermediate_freq','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/meanTb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_std_Tb','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/THot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTHot','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/noise_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/std_dev_noise_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/calibration_time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_sky_elevation_angle','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/TRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TWindow','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdTRoom','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/TOut','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/noise_level','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/VGunn','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/stdVGunn','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/spectrometer1/number_of_hot_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_cold_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)
nccreate(filename,'/spectrometer1/number_of_sky_spectra','Dimensions',{'time',Inf},'Datatype','int64','FillValue',-9999)

nccreate(filename,'/spectrometer1/mean_cold_counts','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_hot_counts','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/spectrometer1/mean_sky_counts','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

%if calibrationTool.savePlanckIntensity
    %nccreate(filename,'/spectrometer1/intensity_planck','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'time',Inf},'Datatype','double','FillValue',-9999);
%end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags dataset

% Coordinates time and flags
nccreate(filename,'/flags/time','Dimensions',{'time',Inf},'Datatype','double')

% if the error vector does not exist, we replace it with a scalar NaN
if isfield(calibratedSpectra,'errorVector')
    lenErrorVect=length(calibratedSpectra(1).errorVector);
else
    lenErrorVect=1;
end
nccreate(filename,'/flags/flags','Dimensions',{'flags',lenErrorVect},'Datatype','int64')

% We input a (1xerrorVectorSize) int vector to identify the errors
nccreate(filename,'/flags/calibration_flags','Dimensions',{'flags',lenErrorVect,'time',Inf},'Datatype','int64','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tipping curve dataset

% Coordinates time and flags
nccreate(filename,'/tipping_curve/time','Dimensions',{'time',Inf},'Datatype','double')

% if the vector does not exist, we replace it with a scalar NaN
if isfield(logFile,'TC')  & calibrationTool.doTippingCurve
    lenTipping=length(logFile.TC(1).sky);
    lenFreqTC=length(logFile.TC(1).tippingCurveRawCounts);
else
    lenTipping=1;
    lenFreqTC = 1;
end
nccreate(filename,'/tipping_curve/tipping_angle','Dimensions',{'tipping_angle',lenTipping,'time',Inf},'Datatype','double')
nccreate(filename,'/tipping_curve/channel_idx_tc','Dimensions',{'channel_idx_tc',lenFreqTC},'Datatype','int64','FillValue',-9999)

nccreate(filename,'/tipping_curve/frequency_tc','Dimensions',{'channel_idx_tc',lenFreqTC},'Datatype','double')

% We input a (1xtipping_angle) int vector to save the tipping curve
% variables
% nccreate(filename,'/tipping_curve/cold_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
% nccreate(filename,'/tipping_curve/hot_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
% nccreate(filename,'/tipping_curve/THot_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/tipping_curve/cold_calib','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/tipping_curve/hot_calib','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/tipping_curve/THot_calib','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/tipping_curve/mean_frequency_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/tipping_curve/opacity_iter_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/tipping_curve/opacity_calib_tc','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
%nccreate(filename,'/tipping_curve/Tb_tc','Dimensions',{'tipping_angle',lenTipping,'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/tipping_curve/cold_spectra','Dimensions',{'channel_idx_tc',lenFreqTC,'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/tipping_curve/cold_spectra_tc','Dimensions',{'tipping_angle',lenTipping,'channel_idx_tc',lenFreqTC,'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/tipping_curve/hot_spectra','Dimensions',{'channel_idx_tc',lenFreqTC,'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/tipping_curve/hot_spectra_tc','Dimensions',{'tipping_angle',lenTipping,'channel_idx_tc',lenFreqTC,'time',Inf},'Datatype','double','FillValue',-9999)

nccreate(filename,'/tipping_curve/sky_spectra','Dimensions',{'tipping_angle',lenTipping,'channel_idx_tc',lenFreqTC,'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meteo Data
nccreate(filename,'/meteo/time','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/air_pressure','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/air_temperature','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/relative_humidity','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)
nccreate(filename,'/meteo/precipitation','Dimensions',{'time',Inf},'Datatype','double','FillValue',-9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write netCDF variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scientific Dataset (spectrometer1,2,...)

% Coordinate variables, directly with their attributes
ncwrite(filename,'/spectrometer1/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/spectrometer1/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/spectrometer1/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/spectrometer1/time','timezone',calibrationTool.timeZone);
ncwriteatt(filename,'/spectrometer1/time','description','mean time recorded at the beginning of all sky measurements during this calibration cycle');

if isfield(calibratedSpectra, 'meanDatetimeMJD2K')
    ncwrite(filename,'/spectrometer1/MJD2K',[calibratedSpectra.meanDatetimeMJD2K]);
    ncwriteatt(filename,'/spectrometer1/MJD2K','units','MJD2K');
    ncwriteatt(filename,'/spectrometer1/MJD2K','calendar','Julian');
    ncwriteatt(filename,'/spectrometer1/MJD2K','description','mean time recorded at the beginning of all sky measurements during this calibration cycle');
end
ncwrite(filename,'/spectrometer1/channel_idx',1:calibrationTool.numberOfChannels);
ncwriteatt(filename,'/spectrometer1/channel_idx','description','index of the spectrometer channels, from 1 to N (number of channels)');

%%%%%%%%%%%%%%%%%
% Geolocation variables
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

%%%%%%%%%%%%%%%%%
% Calibration variables
ncwrite(filename,'/spectrometer1/Tb',vertcat(calibratedSpectra.Tb)');
if calibrationTool.savePlanckIntensity
    %ncwrite(filename,'/spectrometer1/intensity_planck',vertcat(calibratedSpectra.intensityPlanck)');
end  
ncwrite(filename,'/spectrometer1/stdTb',vertcat(calibratedSpectra.stdTb)');
ncwrite(filename,'/spectrometer1/mean_std_Tb',[calibratedSpectra.meanStdTb]);
ncwrite(filename,'/spectrometer1/meanTb',[calibratedSpectra.meanTb]);
ncwrite(filename,'/spectrometer1/frequencies',calibratedSpectra(1).freq);

ncwrite(filename,'/spectrometer1/THot',[calibratedSpectra.THot]);
ncwrite(filename,'/spectrometer1/noise_temperature',[calibratedSpectra.TNoise]);
ncwrite(filename,'/spectrometer1/std_dev_noise_temperature',[calibratedSpectra.stdTNoise]);
ncwrite(filename,'/spectrometer1/calibration_time',60*[calibratedSpectra.calibrationTime]);
ncwrite(filename,'/spectrometer1/mean_sky_elevation_angle',[calibratedSpectra.meanAngleAntenna]);

if isfield(calibratedSpectra,'if')
    ncwrite(filename,'/spectrometer1/intermediate_freq',calibratedSpectra(1).if);
else
    ncwrite(filename,'/spectrometer1/intermediate_freq',vertcat(-9999*ones(length(calibratedSpectra(1).freq),1))');
end

if isfield(calibratedSpectra,'stdTHot')
    ncwrite(filename,'/spectrometer1/stdTHot',[calibratedSpectra.stdTHot]);
else
    ncwrite(filename,'/spectrometer1/stdTHot',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'TempRoom')
    ncwrite(filename,'/spectrometer1/TRoom',[calibratedSpectra.TempRoom]);
    ncwrite(filename,'/spectrometer1/stdTRoom',[calibratedSpectra.stdTempRoom]);
else
    ncwrite(filename,'/spectrometer1/TRoom',-9999*ones(length(calibratedSpectra),1));
    ncwrite(filename,'/spectrometer1/stdTRoom',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'TempOut')
    ncwrite(filename,'/spectrometer1/TOut',[calibratedSpectra.TempOut]);
else
    ncwrite(filename,'/spectrometer1/TOut',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'TempWindow')
    ncwrite(filename,'/spectrometer1/TWindow',[calibratedSpectra.TempWindow]);
else
    ncwrite(filename,'/spectrometer1/TWindow',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'noiseLevel')
    ncwrite(filename,'/spectrometer1/noise_level',[calibratedSpectra.noiseLevel]);
else
    ncwrite(filename,'/spectrometer1/noise_level',-9999*ones(length(calibratedSpectra),1));
end

if isfield(calibratedSpectra,'VGunn')
    ncwrite(filename,'/spectrometer1/VGunn',[calibratedSpectra.VGunn]);
    ncwrite(filename,'/spectrometer1/stdVGunn',[calibratedSpectra.stdVGunn]);
else
    ncwrite(filename,'/spectrometer1/VGunn',-9999*ones(length(calibratedSpectra),1));
    ncwrite(filename,'/spectrometer1/stdVGunn',-9999*ones(length(calibratedSpectra),1));
end

numInd=vertcat(calibratedSpectra.numberOfIndices);
ncwrite(filename,'/spectrometer1/number_of_hot_spectra',numInd(:,1));
ncwrite(filename,'/spectrometer1/number_of_cold_spectra',numInd(:,2));
ncwrite(filename,'/spectrometer1/number_of_sky_spectra',numInd(:,3));

ncwrite(filename,'/spectrometer1/mean_cold_counts',[calibratedSpectra.meanColdCounts]);
ncwrite(filename,'/spectrometer1/mean_hot_counts',[calibratedSpectra.meanHotCounts]);
ncwrite(filename,'/spectrometer1/mean_sky_counts',[calibratedSpectra.meanSkyCounts]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the flags variables 
ncwrite(filename,'/flags/time',[calibratedSpectra.meanDatetime]);
ncwriteatt(filename,'/flags/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/flags/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/flags/time','description','mean time of the measurements for this cycle');
ncwriteatt(filename,'/flags/time','timezone',calibrationTool.timeZone);
if isfield(calibratedSpectra,'errorVector')
    ncwrite(filename,'/flags/flags',1:length(calibratedSpectra(1).errorVector));
    ncwrite(filename,'/flags/calibration_flags',vertcat(calibratedSpectra.errorVector)');
else 
    ncwrite(filename,'/flags/flags',1);
    ncwrite(filename,'/flags/calibration_flags',-9999);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the tipping curve variables 
ncwriteatt(filename,'/tipping_curve/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/tipping_curve/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/tipping_curve/time','description','mean time of the sky observation for this tc');
ncwriteatt(filename,'/tipping_curve/time','timezone',calibrationTool.timeZone);
if isfield(logFile,'TC') & calibrationTool.doTippingCurve
    ncwrite(filename,'/tipping_curve/time',[logFile.TC.meanDateNum]);
    ncwrite(filename,'/tipping_curve/channel_idx_tc',logFile.TC(1).channels);
    ncwrite(filename,'/tipping_curve/tipping_angle',vertcat(logFile.TC.skyAngle)');
    %ncwrite(filename,'/tipping_curve/sky',vertcat(logFile.TC.sky)');
    %ncwrite(filename,'/tipping_curve/Tb_tc',vertcat(logFile.TC.Tb_Calib)');
    ncwrite(filename,'/tipping_curve/hot_spectra',vertcat(logFile.TC.hotSpectra)');
    ncwrite(filename,'/tipping_curve/cold_spectra',vertcat(logFile.TC.coldSpectra)');
    ncwrite(filename,'/tipping_curve/hot_calib', [logFile.TC.hotCalib]);
    ncwrite(filename,'/tipping_curve/THot_calib', [logFile.TC.THotCalib]);
    ncwrite(filename,'/tipping_curve/cold_calib', [logFile.TC.coldCalib]);
    ncwrite(filename,'/tipping_curve/hot_calib', [logFile.TC.hotCalib]);
    ncwrite(filename,'/tipping_curve/opacity_calib_tc', [logFile.TC.tauCalibZenith]);
    if isfield(logFile.TC,'')
        ncwrite(filename,'/tipping_curve/opacity_iter_tc', [logFile.TC.tauIter]);
    end
    ncwrite(filename,'/tipping_curve/mean_frequency_tc', [logFile.TC.meanFreq]);
    ncwrite(filename,'/tipping_curve/frequency_tc', logFile.TC(1).frequency);   
    for j=1:length(logFile.TC)
        ncwrite(filename,'/tipping_curve/sky_spectra',logFile.TC(j).sky_spectra,[1,1,j]);
        if isfield(logFile.TC,'cold_spectraTC')
            ncwrite(filename,'/tipping_curve/hot_spectra_tc',logFile.TC(j).hot_spectraTC,[1,1,j]);
            ncwrite(filename,'/tipping_curve/cold_spectra_tc',logFile.TC(j).cold_spectraTC,[1,1,j]);
        end
    end
else
    ncwrite(filename,'/tipping_curve/time',calibratedSpectra(1).firstSkyTime);
    ncwrite(filename,'/tipping_curve/channel_idx_tc',-9999);
    ncwrite(filename,'/tipping_curve/tipping_angle',-9999);
    %ncwrite(filename,'/tipping_curve/sky',-9999);
    %ncwrite(filename,'/tipping_curve/Tb_tc',-9999);
    ncwrite(filename,'/tipping_curve/hot_calib', -9999);
    ncwrite(filename,'/tipping_curve/THot_calib', -9999);
    ncwrite(filename,'/tipping_curve/cold_calib', -9999);
    ncwrite(filename,'/tipping_curve/hot_calib', -9999);
    ncwrite(filename,'/tipping_curve/opacity_calib_tc', -9999);
    ncwrite(filename,'/tipping_curve/mean_frequency_tc', -9999);
    ncwrite(filename,'/tipping_curve/frequency_tc', -9999);
    ncwrite(filename,'/tipping_curve/opacity_iter_tc', -9999)
%     ncwrite(filename,'/tipping_curve/cold_tc', -9999);
%     ncwrite(filename,'/tipping_curve/hot_tc', -9999);
%     ncwrite(filename,'/tipping_curve/hot_tc', -9999);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing Meteo Data
if isfield(logFile,'meteo') && ~isempty(fieldnames(logFile.meteo))
    ncwrite(filename,'/meteo/time',[logFile.meteo.dateNum]);
    ncwrite(filename,'/meteo/air_pressure',[logFile.meteo.air_pressure]);
    ncwrite(filename,'/meteo/air_temperature',[logFile.meteo.air_temperature]);
    ncwrite(filename,'/meteo/relative_humidity',[logFile.meteo.rel_humidity]);
    ncwrite(filename,'/meteo/precipitation',[logFile.meteo.precipitation]);
else
    ncwrite(filename,'/meteo/time',-9999);
    ncwrite(filename,'/meteo/air_pressure',-9999);
    ncwrite(filename,'/meteo/air_temperature',-9999);
    ncwrite(filename,'/meteo/relative_humidity',-9999);
    ncwrite(filename,'/meteo/precipitation',-9999);
end

% add time attributes directly
ncwriteatt(filename,'/meteo/time','units',calibrationTool.meanDatetimeUnit);
ncwriteatt(filename,'/meteo/time','calendar',calibrationTool.calendar);
ncwriteatt(filename,'/meteo/time','description','time from the meteo stations');
ncwriteatt(filename,'/meteo/time','timezone',calibrationTool.timeZone);

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
    %ncwriteatt(filename,'/flags/calibration_flags','errorCode_1',calibratedSpectra(1).errorVectorDescription(1));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_2',calibratedSpectra(1).errorVectorDescription(2));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_3',calibratedSpectra(1).errorVectorDescription(3));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_4',calibratedSpectra(1).errorVectorDescription(4));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_5',calibratedSpectra(1).errorVectorDescription(5));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_6',calibratedSpectra(1).errorVectorDescription(6));
%     ncwriteatt(filename,'/flags/calibration_flags','errorCode_7',calibratedSpectra(1).errorVectorDescription(7));
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

%attrVal.effCalTime = {'effective calibration time',...
%    'calibration_time',...
%    'second',...
 %   'still to improve'};

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% attrs or TC group
attrVal.tipping_angle = {'tipping curve elevation angles',...
    'tipping_curve_angle',...
    'degree',...
    'elevation angles used for the tipping curve'};

attrVal.cold_calib = {'mean FFTS cold counts for tipping curve',...
    'mean_cold_counts_tc',...
    '1',...
    'mean raw FFTS counts on cold load during this cycle'};

attrVal.hot_calib = {'mean FFTS hot counts for tipping curve',...
    'mean_hot_counts_tc',...
    '1',...
    'mean raw FFTS counts on hot load during this cycle'};

attrVal.sky = {'mean FFTS sky counts for tipping curve',...
    'mean_sky_counts_tc',...
    '1',...
    'mean raw FFTS counts on sky during TC'};

attrVal.sky_spectra = {'FFTS sky counts for tipping curve',...
    'sky_counts_tc',...
    '1',...
    'raw FFTS counts on sky during TC'};

attrVal.tau_tc = {'tau TC',...
    'tau_tc',...
    '1',...
    'estimated tau from TC'};

attrVal.THot_calib = {'hot temperature calibration',...
    'THot_calib',...
    'K',...
    'mean temperature of the hot load during calibration cycle around TC'};

attrVal.hot_spectra = {'hot spectra calibration',...
    'hot_spectra',...
    '1',...
    'raw FFTS counts on hot during calibration cycle around TC'};

attrVal.hot_spectra_tc = {'hot spectra TC',...
    'hot_spectra_tc',...
    '1',...
    'raw FFTS counts on hot during TC'};

attrVal.cold_spectra = {'cold spectra calibration',...
    'cold_spectra',...
    '1',...
    'raw FFTS counts on cold during calibration cycle around TC'};

attrVal.cold_spectra_tc = {'cold spectra TC',...
    'cold_spectra_tc',...
    '1',...
    'raw FFTS counts on cold during TC'};

attrVal.frequency_tc = {'frequency tc',...
    'frequency_tc',...
    'Hz',...
    'frequency vector used for the tc calibration'};

attrVal.mean_frequency_tc = {'mean frequency tc',...
    'frequency_tc',...
    'Hz',...
    'mean frequency used for the tc calibration'};

attrVal.opacity_iter_tc = {'opacity iter tc',...
    'opacity_iter_tc',...
    '-',...
    'opacity at zenith iterated with tc measurements only'};

attrVal.opacity_calib_tc = {'mean frequency tc',...
    'opacity_calib_tc',...
    '-',...
    'opacity at zenith computed using the calibration cycle around tc'};

% Ugly and open to suggestion
for i=1:length(attrName)
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
    ncwriteatt(filename,'/spectrometer1/frequencies',attrName{i},attrVal.freq{i});
    ncwriteatt(filename,'/spectrometer1/intermediate_freq',attrName{i},attrVal.if{i});
    ncwriteatt(filename,'/spectrometer1/THot',attrName{i},attrVal.THot{i});
    ncwriteatt(filename,'/spectrometer1/stdTHot',attrName{i},attrVal.stdTHot{i});
    ncwriteatt(filename,'/spectrometer1/noise_temperature',attrName{i},attrVal.TNoise{i});
    ncwriteatt(filename,'/spectrometer1/std_dev_noise_temperature',attrName{i},attrVal.stdTNoise{i});
    ncwriteatt(filename,'/spectrometer1/calibration_time',attrName{i},attrVal.calibrationTime{i});
    ncwriteatt(filename,'/spectrometer1/mean_sky_elevation_angle',attrName{i},attrVal.meanAngleAntenna{i});
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


    
    % Meteo attr
    ncwriteatt(filename,'/meteo/air_pressure',attrName{i},attrVal.airP{i});
    ncwriteatt(filename,'/meteo/air_temperature',attrName{i},attrVal.airT{i});
    ncwriteatt(filename,'/meteo/relative_humidity',attrName{i},attrVal.relH{i});
    ncwriteatt(filename,'/meteo/precipitation',attrName{i},attrVal.precipitation{i});
    
    % TC attrs
    ncwriteatt(filename,'/tipping_curve/tipping_angle',attrName{i},attrVal.tipping_angle{i});
    ncwriteatt(filename,'/tipping_curve/cold_calib',attrName{i},attrVal.cold_calib{i});
    ncwriteatt(filename,'/tipping_curve/hot_calib',attrName{i},attrVal.hot_calib{i});
    ncwriteatt(filename,'/tipping_curve/THot_calib',attrName{i},attrVal.THot_calib{i});
    ncwriteatt(filename,'/tipping_curve/sky_spectra',attrName{i},attrVal.sky_spectra{i});
    ncwriteatt(filename,'/tipping_curve/hot_spectra',attrName{i},attrVal.hot_spectra{i});
    ncwriteatt(filename,'/tipping_curve/cold_spectra',attrName{i},attrVal.cold_spectra{i});
    ncwriteatt(filename,'/tipping_curve/hot_spectra_tc',attrName{i},attrVal.hot_spectra_tc{i});
    ncwriteatt(filename,'/tipping_curve/cold_spectra_tc',attrName{i},attrVal.cold_spectra_tc{i});
    ncwriteatt(filename,'/tipping_curve/opacity_calib_tc',attrName{i},attrVal.opacity_calib_tc{i});
    ncwriteatt(filename,'/tipping_curve/opacity_iter_tc',attrName{i},attrVal.opacity_iter_tc{i});
    ncwriteatt(filename,'/tipping_curve/mean_frequency_tc',attrName{i},attrVal.mean_frequency_tc{i});
    ncwriteatt(filename,'/tipping_curve/frequency_tc',attrName{i},attrVal.frequency_tc{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding debug group (optionnal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are 2 types of debug information for which we create new groups.
if calibrationTool.calType=="debug"
    nccreate(filename,'/UpDown/timeUpDown','Dimensions',{'timeUpDown',Inf},'Datatype','double','Format','netcdf4');
    nccreate(filename,'/UpDown/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','int64','FillValue',-9999)
    
    nccreate(filename,'/UpDown/TbUp','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'timeUpDown',Inf},'Datatype','double','FillValue',-9999);
    nccreate(filename,'/UpDown/TbDown','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'timeUpDown',Inf},'Datatype','double','FillValue',-9999);

    ncwrite(filename,'/UpDown/TbUp',vertcat(calibratedSpectra.TbUp)');
    ncwrite(filename,'/UpDown/TbDown',vertcat(calibratedSpectra.TbDown)');
    
    % For all cycle:
    % initialize matrices
    nClean=0;
    for t = 1:length(calibratedSpectra)
        nClean=nClean+size(calibratedSpectra(t).TbAll,1);
    end
    
    channel_idx=int64(ones(nClean,calibrationTool.numberOfChannels)*NaN);
    TbAll=single(ones(nClean,calibrationTool.numberOfChannels)*NaN);
    cycleNumber=int64(ones(nClean,1)*NaN);
    cycleId=int64(1:nClean);
    
    counter=1;
    for t = 1:length(calibratedSpectra)
        for i = 1:size(calibratedSpectra(t).TbAll,1)
            TbAll(counter,:)=calibratedSpectra(t).TbAll(i,:);
            channel_idx(counter,:)=1:calibrationTool.numberOfChannels;
            cycleNumber(counter)=t;
            counter=counter+1;
        end
    end
    
    % Group for debugging variables:
    nccreate(filename,'/debug/cycle_id','Dimensions',{'cycle_id',Inf},'Datatype','int64')
    nccreate(filename,'/debug/channel_idx','Dimensions',{'channel_idx',calibrationTool.numberOfChannels},'Datatype','int64','FillValue',-9999)
    
    nccreate(filename,'/debug/calibration_cycle_number','Dimensions',{'cycle_id',Inf},'Datatype','int64','FillValue',-9999)
    nccreate(filename,'/debug/Tb_all','Dimensions',{'channel_idx',calibrationTool.numberOfChannels,'cycle_id',Inf},'Datatype','single','FillValue',-9999)
    
    % Writing the debug variables
    ncwrite(filename,'/debug/cycle_id',cycleId);
    
    ncwrite(filename,'/debug/channel_idx',1:calibrationTool.numberOfChannels);
    ncwrite(filename,'/debug/calibration_cycle_number',cycleNumber');
    ncwrite(filename,'/debug/Tb_all',TbAll');   
end

disp(['File saved as: ' filename])
end
