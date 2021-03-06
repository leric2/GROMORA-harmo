function meteoData = read_meteo_data_payerne(calibrationTool)
%==========================================================================
% NAME          | read_meteo_data_payerne
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez & Eliane Maillard (adapted to GROSOM by ES)
% CREATION      | 09.2014 // modified 03.2020
%               |
% ABSTRACT      | Function to read temperature and humidity values measured
%               | by the ASTA station from MCH in Payerne. It is then saved
%               | in the form of a standard meteoData structure (see
%               | outputs). If the local meteo station file does not exist,
%               | we try reading the anetz meteo data directly. 
%               | 
%               |
% ARGUMENTS     | INPUTS:   1. calibrationTool:
%               |               - dateStr
%               |               - meteoFolder
%               |               - referenceTime
%               |               - timeZone
%               |               - zeroDegInKelvin
%               |               - meteoAnetzFolder
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
    % fileName for meteo file in Payerne: yyyyjjj.002
    % julian day
    jjj= datenum(calibrationTool.dateStr)-datenum(str2num(datestr(datenum(calibrationTool.dateStr),'yyyy')),01,01);
    
    % First reading the Meteo dataset for this day
    MetNum=str2num(calibrationTool.dateStr(1:4))*1000;
    dateStringMeteo=[num2str(MetNum+jjj+1) '.002'];
    %dateStringMeteo='2020065.002';
    
    meteoDataFile=[calibrationTool.meteoFolder dateStringMeteo];
    
    % Transforming it into matlab structure
    fid=fopen(meteoDataFile,'r');
    met=fscanf(fid,'%f %f %f %f %f %f',[6 inf]);
    met=met';
    
    met(met == -999.0) = NaN;
    met(met == 3276.7) = NaN;
    
    datemet=datevec(met(:,1)-1+datenum(str2num(dateStringMeteo(1:4)),01,01));
    
    dateNum=datenum(datemet)-calibrationTool.referenceTime;
    dateTime=datetime(datemet);
    
    meteoData=struct();
    
    for i = 1:length(met)
        meteoData(i).dateNum=dateNum(i);
        meteoData(i).dateTime=dateTime(i);
        meteoData(i).dateTime.TimeZone =calibrationTool.timeZone;
        meteoData(i).air_temperature=met(i,2) + calibrationTool.zeroDegInKelvin;
        meteoData(i).rel_humidity=met(i,3);
        meteoData(i).air_pressure=met(i,4);
        meteoData(i).precipitation=met(i,5);
        %meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
        meteoData(i).tod = 24*(met(i,1) - met(1,1));
    end
catch ME
    warning('no meteo from local station, trying with anetz data')
    try
        meteoData = read_meteo_data_MCH(calibrationTool);
        %disp('Found anetz data for this day');
    catch ME
        warning(ME.message)
        disp('no meteo data loaded for this day')
        meteoData = struct();
    end
end
