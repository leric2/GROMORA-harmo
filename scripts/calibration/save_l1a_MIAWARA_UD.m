function calibrationTool = save_l1a_MIAWARA_UD(calibrationTool,logFile,calibratedSpectra,warningLevel0)
%==========================================================================
% NAME      | save_level1a_daily.m
% TYPE      | function
% AUTHOR(S) | Alistair Bell
% CREATION  | 03.2023
% ABSTRACT  | Updated version of the write netCDF function
%           | originally written by Eric Sauvageat. Instead
%           | of writing all text inside this function, the 
%           | text is imported from a csv file. The function
%           | tests each key in the csv file to see if 
%           | it exists in the specified structure function. If 
%           | it is contained, the field will be created, the 
%           | data written to the netCDF file and the attributes
%           | which are also contained in the netCDF file will
%           | be appended.
%           | Accepted dimension sizes currently are 'ChanNos' or 'Inf',
%           | any other dimensions will cause a problem unless this file 
%           | is correctly edited. 
%           |
%           |
%           |
%           |The function is also adapted to multiple spectrometer- if
%           |present, this must be specified in CalibrationTool, along with
%           |the number of channels that each spectrometer has, which
%           |should be detailed as an array in calibrationTool.spectroSubChan
%           |where each element corresponds to that spectrometer.
%           |
%           |
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

%Temporary File name TO DELETE
%filename = '/home/alistair/Downloads/work/tempdata/temp_calib.nc';
filename = calibrationTool.filenameLevel1a;

%filename of lookup table containing all relevent netcdf info
filename_lookup = calibrationTool.filenameLookup;
filename_meteo_lookup = calibrationTool.filenameMeteoLookup;

% Rewrite if already existing
if  exist(filename,'file')
    delete(filename)
end

%lookup table for netCDF fields
opts = detectImportOptions(filename_lookup);
opts.VariableTypes = {'string', 'string', 'string','string','string','string','double', 'string', 'string', 'string', 'string'};
lookUpFields = readtable('/home/alistair/MIAWARA_ret/extra_files/netCDF_fields_LN2.csv', opts);

opts = detectImportOptions(filename_meteo_lookup);
opts.VariableTypes = {'string', 'string', 'string','string','double', 'string', 'string', 'string', 'string'};
lookUpFieldsMeteo = readtable(filename_meteo_lookup, opts);

%make the dimensions into a cell array
temdims = lookUpFields.Dimensions;
csvContent = erase(temdims, {'{', '}'});
rows = split(csvContent, ';');
numRows = numel(temdims);

% Initialize the cell array
cellArray = cell(numRows, 1);

%iterate through rows and turn into cell
for i = 1:length(cellArray)
    rowContent = split(rows(i), ',');
    cellArray(i) = {rowContent'};
end

%assign new cell array as dimensions
lookUpFields.Dimensions = cellArray;

%if sepectrometerQuantity not set, assume 1
if ~isfield(calibrationTool, 'spectrometerQuantity')
    calibrationTool.spectrometerQuantity = 1;
    calibrationTool.spectroSubChan = [calibrationTool.numberOfChannels];
end

initChan = 1;
for j = 1:calibrationTool.spectrometerQuantity
    stChan(j) = initChan;
    endChan(j)= stChan(j)+calibrationTool.spectroSubChan(j)-1;
    initChan = endChan(j)+1;
end

%define the channel indices 
%could be improved i.e. added as an earlier variable in calibrationTool
calibrationTool.chanIdx = [1:calibrationTool.numberOfChannels];

%all saved arrays must be single column - need to adapt frequency field
for i = 1:length(calibratedSpectra)
    calibratedSpectra(i).frequency = calibratedSpectra(i).freq(:)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scientific Dataset (spectrometer1,2,...)

attrName = {'long_name','standard_name','units','description'};

for i = 1:calibrationTool.spectrometerQuantity
    % Create dynamic field name
    N_spectro = ['spectrometer' num2str(i)];
    chanNos = calibrationTool.spectroSubChan(i);
    
  
   for k = 1:height(lookUpFields)
       %select the structure from which data will be exported
       tempStruct = eval(lookUpFields.StructName{k});

       %check that field exists in this structure
       if isfield(tempStruct, lookUpFields.CS_name{k})

           %select dimensions of the field
           tempDims = cellstr(lookUpFields.Dimensions{k,1});

           for l = 2:2: length(tempDims)
             %for each dimension in the csv file, assign numerical value 
             switch tempDims{l}
               case 'Inf'
                 disp('gottit inf')
                 tempDims{l} = Inf;
               case 'chanNos'
                 disp('gottit chanNos')
                 tempDims{l} = chanNos;
               end
           end
          
           %lookUpFields.Dimensions(k) = {tempDims};
           ncFieldName = strcat( N_spectro, '/' , lookUpFields.NC_name(k));
           nccreate(filename,ncFieldName,'Dimensions', tempDims,'Datatype',lookUpFields.Datatype(k),'FillValue', lookUpFields.FillValue(k), 'Format', lookUpFields.Format(k));
           numDims = length(tempDims)/2;

           %get array from structure
           tempArray = vertcat(tempStruct.(lookUpFields.CS_name(k)));

           switch numDims
               case 0
                   %add something here
               case 1
                 if strcmp(lookUpFields.Dimensions{k}{2}, "Inf")
                    if length(tempArray) == length(calibratedSpectra)
                        ncwrite(filename,ncFieldName,tempArray)
                    elseif length(tempArray) == 1
                        tempArray  = tempArray * ones(length(calibratedSpectra),1);
                        ncwrite(filename,ncFieldName,tempArray)
                    else
                        fprintf('Inconsistant Dimensions for variable: %s \n', lookUpFields.CS_name(k));
                    end

                 elseif strcmp(lookUpFields.Dimensions{k}{2}, 'chanNos')
                     if length(tempArray) >= endChan
                         ncwrite(filename,ncFieldName,tempArray(stChan:endChan))
                     elseif length(tempArray) == 1
                         ncwrite(filename,ncFieldName,tempArray*ones(chanNos,1))
                     else
                        fprintf('Inconsistant Dimensions for variable: %s \n', lookUpFields.CS_name(k));
                     end
                 end
               case 2
                 if length(size(tempArray))==2 & any(size(tempArray)>1)
                     if strcmp(lookUpFields.Dimensions{k}{2}, Inf)
                       ncwrite(filename,ncFieldName,tempArray(:,stChan(i):endChan(i)))
                     else
                       ncwrite(filename,ncFieldName,tempArray(:, stChan(i):endChan(i) )' )
                     end
                 else
                   fprintf('Inconsistant Dimensions for variable: %s \n', lookUpFields.CS_name(k));
                 end

              case 3
                 %currently should be nothing to append with three dimensions
           end
           for j = 1:length(attrName)
             ncwriteatt(filename,ncFieldName,attrName{j},lookUpFields.(attrName{j}){k});
           end
      
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Meteo Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing Meteo Data
if isfield(logFile,'meteo') && ~isempty(fieldnames(logFile.meteo))
    %writetime = convertTo(writetime, 'epochtime','Epoch','1970-01-01')/(24.*60.*60.);
    %ncwrite(filename,'/meteo/time',writetime);
    disp(length(logFile.meteo));
    for i = 1:length(logFile.meteo)
        logFile.meteo(i).time_sSince1970 = 86400*(logFile.meteo(i).dateNum - datenum(1970,1,1));
    end
    
    for k = 1:height(lookUpFieldsMeteo)
        
     %check that field exists in this structure
      if isfield(logFile.meteo, lookUpFieldsMeteo.LF_meteo_name{k})
        ncFieldName = strcat( '/meteo/' , lookUpFieldsMeteo.NC_name(k));
        nccreate(filename,ncFieldName,'Dimensions',{'time',Inf},'Datatype',lookUpFields.Datatype(k),'FillValue', lookUpFields.FillValue(k), 'Format', lookUpFields.Format(k));
        tempArray = vertcat(logFile.meteo.(lookUpFieldsMeteo.LF_meteo_name{k}));
        ncwrite(filename,ncFieldName,tempArray)
        for j = 1:length(attrName)
             ncwriteatt(filename,ncFieldName,attrName{j},lookUpFieldsMeteo.(attrName{j}){k});
        end
      end
    end
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
ncwriteatt(filename,'/','filename',filename);
ncwriteatt(filename,'/','creation_date',datestr(datetime('now','TimeZone','Z'),'yyyymmddTHHMMSSZ'));
ncwriteatt(filename,'/','featureType','timeSeries');


disp(['File saved as: ' filename])
end
