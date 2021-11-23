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
    % ATTENTION, day must be real day + 5.
    
    % First reading the Meteo dataset for this day
    dateStringMeteo=['VQCA43.' datestr(calibrationTool.dateTime + days(5), 'YYYYmmDD') '0200'];
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
        meteoData(i).air_temperature = meteoData(i).air_temperature + calibrationTool.zeroDegInKelvin;
    end
catch ME
    warning(ME.message)
    try 
        dateStringMeteo=['meteo_' calibrationTool.anetzStnName '.csv'];
    
        meteoDataFile=[calibrationTool.meteoAnetzExtraFolder dateStringMeteo];
    
        % Transforming it into matlab structure
        T=readtable(meteoDataFile,'HeaderLines',1 , 'FileType','text');
        T.Properties.VariableNames = {'year' 'month' 'day' 'hour' 'minute',...
            'second' 'stationid' 'air_pressure' 'rel_humidity' 'air_temperature' 'grass_temperature',...
            'precipitation' 'wind_dirn'  'wind_speed' 'max_gust' 'distant_lightning' 'near_lightning' 'global_rad'};
        dateN = datenum(T.year, T.month, T.day, T.hour, T.minute, T.second);

        % For speed reason, we keep only the 2 days before and after
        % (useless anyway)
        T.dateNum = dateN;
        T = T(T.dateNum > calibrationTool.timeNumber-2 & T.dateNum < calibrationTool.timeNumber+2, :);
        meteoData=table2struct(T);
        
        
        for i = 1:length(meteoData)
            meteoData(i).dateNum=T.dateNum(i)-calibrationTool.referenceTime;
            meteoData(i).dateTime=datetime(T.dateNum(i),'ConvertFrom','datenum');
            meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
            %meteoData(i).dateTime=datetime(meteoData(i).time,'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSSSSSZ');
            meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
            meteoData(i).air_temperature = meteoData(i).air_temperature + calibrationTool.zeroDegInKelvin;
        end
    catch ME
        disp('no MCH data loaded for this day')
        meteoData = struct();
end
end
