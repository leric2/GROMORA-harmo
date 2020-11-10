function [data, meteoData, calibrationTool] = read_level1b_daily(calibrationTool)
%==========================================================================
% NAME          | read_level1a_daily.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Function to read a level1a previously saved.
%               |
%               |
%               |
% ARGUMENTS     | INPUTS: -retrievalTool (containing the filename of the
%               | level1a file to read.
%               |
%               | OUTPUTS: - calib (part of "calibratedSpectra")
%               |          - retrievalTool
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================

% filename:
filename=calibrationTool.filenameLevel1b;

ncinfo(filename)
gNames = {ncinfo(filename).Groups.Name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = struct();
meteoData= struct();
for g=1:length(gNames)
    gName = gNames{g};
    
    if strcmp(gName,'meteo')
        vNames = {ncinfo(filename,gName).Variables.Name};
        %dimNames = {ncinfo(filename,gName).Dimensions.Name};
        %         size=zeros(1,length(dimNames));
        %         for i=1:length(dimNames)
        %             size(i) = ncinfo(filename,fullfile(gName,dimNames{i})).Size;
        %         end
        
        for v=1:length(vNames)
            varName = vNames{v};
            
            if length(ncinfo(filename,fullfile(gName,varName)).Dimensions)==1
                if strcmp(ncinfo(filename,fullfile(gName,varName)).Dimensions.Name,'time')
                    meteoData.(varName) = ncread(filename,fullfile(gName,varName))';
                else
                    meteoData(t).(varName) = ncread(filename,fullfile(gName,varName),1,Inf);
                end
            else
                meteoData(t).(varName)= ncread(filename,fullfile(gName,varName),[1,t],[Inf,1])';
            end
        end
    else
        vNames = {ncinfo(filename,gName).Variables.Name};
        dimNames = {ncinfo(filename,gName).Dimensions.Name};
        size=zeros(1,length(dimNames));
        for i=1:length(dimNames)
            size(i) = ncinfo(filename,fullfile(gName,dimNames{i})).Size;
        end
        
        % Intialize our structure with time
        for t = 1:size(1)
            varName = vNames{1};
            data(t).(varName) = ncread(filename,fullfile(gName,varName),t,1);
        end
        dateT = num2cell(datetime([data.time] + datenum(1970,1,1),'ConvertFrom','datenum'));
        [data.meanDateTime] = dateT{:};
        %if strcmp('time',dimName)
        for v=2:length(vNames)
            varName = vNames{v};
            
            if length(ncinfo(filename,fullfile(gName,varName)).Dimensions)==1
                if strcmp(ncinfo(filename,fullfile(gName,varName)).Dimensions.Name,'time')
                    dat=num2cell(ncread(filename,fullfile(gName,varName)));
                    [data.(varName)] = dat{:};
                else
                    for t = 1:size(1)
                        data(t).(varName) = ncread(filename,fullfile(gName,varName),1,Inf)';
                    end
                end
            else
                %dat=num2cell(ncread(filename,fullfile(gName,varName)))';
                for t = 1:size(1)
                    data(t).(varName)= ncread(filename,fullfile(gName,varName),[1,t],[Inf,1])';
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading attributes
calibrationTool.logFile.raw_file_warning=ncreadatt(filename,'/','raw_file_warning');
calibrationTool.logFile.comment=ncreadatt(filename,'/','comment');
calibrationTool.logFile.raw_file_comment=ncreadatt(filename,'/','raw_file_comment');
calibrationTool.logFile.rawFilename=ncreadatt(filename,'/','raw_filename');
calibrationTool.logFile.rawData=ncreadatt(filename,'/','raw_data');
calibrationTool.logFile.raw_data_software_version=ncreadatt(filename,'/','raw_data_software_version');
calibrationTool.logFile.calibration_version=ncreadatt(filename,'/','calibration_version');
calibrationTool.logFile.creation_date_level1a=ncreadatt(filename,'/','creation_date');
calibrationTool.logFile.raw_data_software_version=ncreadatt(filename,'/','raw_data_software_version');
calibrationTool.logFile.filenameLevel1a=ncreadatt(filename,'/','filename');

calibrationTool.flagVectorLength = ncreadatt(filename,'/flags/','number_of_flags');

% Coordinate variables, directly adding the attributes
calibrationTool.timeUnit = ncreadatt(filename,'/spectrometer1/time','units');
calibrationTool.timeCalendar= ncreadatt(filename,'/spectrometer1/time','calendar');


meteoData.dateTime = datetime(meteoData.time + datenum(1970,1,1),'ConvertFrom','datenum');

if ~isfield(data,'noise_level')
     warning('no noise level defined for calibration');
     in = num2cell(-9999*ones(1,length(data)));
     [data.noiseLevel] = in{:};
end

disp(['File read : ' filename])

end


