function meteoData = read_meteo_data_MCH(calibrationTool)
%==========================================================================
% NAME          | read_meteo_data_MCH
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 03.2021
%               |
% ABSTRACT      | Function to read temperature and humidity values measured
%               | by the MCH ground stations.
%               | It is then saved
%               | in the form of a standard meteoData structure (see
%               | outputs).
%               | 
%               |
% ARGUMENTS     | INPUTS:   1. calibrationTool:
%               |               - dateStr
%               |               - meteoAnetzFolder
%               |               - referenceTime
%               |               - timeZone
%               |               - anetzStnName
%               |
%               | OUTPUTS:  1. meteoData: standard structure containing:
%               |               - dateTime
%               |               - dateNum
%               |               - air_temperature
%               |               - rel_humidity
%               |               - air_pressure
%               |               - precipitation
%               |               - tod
%               |
%==========================================================================
%%% CHECK FOR LEAP YEARS ?
try
    % fileName for meteo file in anetz: VQCA43.YYYYMMDD0200
    
    % First reading the Meteo dataset for this day
    dateStringMeteo=['VQCA43.' calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10) '0200'];
    %dateStringMeteo='2020065.002';
    
    meteoDataFile=[calibrationTool.meteoAnetzFolder dateStringMeteo];
    
    % Transforming it into matlab structure
    T=readtable(meteoDataFile,'HeaderLines',4 , 'FileType','text');
    T.Properties.VariableNames = {'stn' 'time' 'winddirn' 'windspeed' 'maxgust',...
        'air_temperature' 'grasstemp5cm' 'rel_humidity' 'air_pressure' 'precipitation',...
        'nahblitz' 'fernblitz'  'globalrad' 'sunduration' 'skybright_current' 'skybright_avg'};
    local = T(find(strcmp(calibrationTool.anetzStnName,T.stn)),:);
    
    meteoData=table2struct(local);
    
    for i = 1:length(meteoData)
        dateN = datenum(string(meteoData(i).time),'YYYYmmddHHMMZ');
        meteoData(i).dateNum=dateN-calibrationTool.referenceTime;
        meteoData(i).dateTime=datetime(dateN,'ConvertFrom','datenum');
        meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
        %meteoData(i).dateTime=datetime(meteoData(i).time,'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSSSSSZ');
        meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
    end
catch ME
    warning(ME.message)
    disp('no MCH data loaded for this day')
    meteoData = struct();
end
end
