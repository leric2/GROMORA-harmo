function calibratedSpectra = get_meteo_data_payerne(calibratedSpectra,retrievalTool)
%==========================================================================
% NAME          | get_meteo_data_payerne
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez & Eliane Maillard (adapted to GROSOM by ES)
% CREATION      | 09.2014 // modified 03.2020
%               |
% ABSTRACT      | Function to read temperature and humidity values measured
%               | by the ASTA station from MCH in Payerne
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

%%% CHECK FOR LEAP YEARS ?

% fileName for meteo file in Payerne: yyyyjjj.002
% julian day
jjj= datenum(retrievalTool.dateStr)-datenum(str2num(datestr(datenum(retrievalTool.dateStr),'yyyy')),01,01);

% First reading the Meteo dataset for this day
MetNum=str2num(retrievalTool.dateStr(1:4))*1000;
dateStringMeteo=[num2str(MetNum+jjj+1) '.002'];
%dateStringMeteo='2020065.002';

meteoDataFile=[retrievalTool.meteoFolder dateStringMeteo];

%opts = detectImportOptions(meteoDataFile)

% Transforming it into matlab structure
fid=fopen(meteoDataFile,'r');
met=fscanf(fid,'%f %f %f %f %f %f',[6 inf]);
met=met';

met(met == -999.0) = NaN;
met(met == 3276.7) = NaN;

datemet=datevec(met(:,1)-1+datenum(str2num(dateStringMeteo(1:4)),01,01));
dateTime=datenum(datemet)-datenum(1970,1,1);

meteoData=struct();

for i = 1:length(met)
    meteoData(i).dateTime=dateTime(i);
    meteoData(i).air_temperature=met(i,2);
    meteoData(i).rel_humidity=met(i,3);
    meteoData(i).air_pressure=met(i,4);
    meteoData(i).precipitation=met(i,5);
end

for t=1:length(calibratedSpectra)
    %start=datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum');
    stop=datenum(datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum')+seconds(calibratedSpectra(t).calibrationTime))-datenum(1970,1,1);
    
    % Selecting the interesting values for each calibration cycle:
    rowInd=([meteoData.dateTime]>=calibratedSpectra(t).timeMin & [meteoData.dateTime]<=stop);
    
    calibratedSpectra(t).meanAirTemperature=nanmean(vertcat(meteoData(rowInd).air_temperature))+273.15;
    calibratedSpectra(t).meanRelHumidity=nanmean(vertcat(meteoData(rowInd).rel_humidity));
    calibratedSpectra(t).meanAirPressure=nanmean(vertcat(meteoData(rowInd).air_pressure));
    calibratedSpectra(t).rainAccumulation=nansum(vertcat(meteoData(rowInd).precipitation));
end
    
end
