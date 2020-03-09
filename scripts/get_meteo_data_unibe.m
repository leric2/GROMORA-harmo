function calibratedSpectra = get_meteo_data_unibe(calibratedSpectra,retrievalTool)
%==========================================================================
% NAME          | get_meteo_data_unibe
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014 // modified 01.2020
%               |
% ABSTRACT      | Function to read temperature and humidity values measured
%               | by the sensors on the top of ExWi building at Unibe
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

if datetime(calibratedSpectra(1).year,calibratedSpectra(1).month,calibratedSpectra(1).day) > datetime(2017,08,10)
    % First reading the Meteo dataset for this day
    dateStringMeteo=[retrievalTool.dateStr(1:4) '-' retrievalTool.dateStr(6:7) '-' retrievalTool.dateStr(9:10)];
    meteoDataFile=[retrievalTool.meteoFolder 'meteo_' dateStringMeteo '.csv'];
    
    % Transforming it into matlab structure
    T=readtable(meteoDataFile);
    meteoData=table2struct(T);
    
    for i = 1:length(meteoData)
        meteoData(i).dateTime=datenum(meteoData(i).time,'yyyy-mm-ddTHH:MM:SS.FFFZ')-datenum(1970,1,1);
    end
    
    %%%%%% Storing values into calibratedSpectra
    for t=1:length(calibratedSpectra)
        %start=datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum');
        stop=datenum(datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum')+seconds(calibratedSpectra(t).calibrationTime))-datenum(1970,1,1);
        
        % Selecting the interesting values for each calibration cycle:
        rowInd=([meteoData.dateTime]>=calibratedSpectra(t).timeMin & [meteoData.dateTime]<=stop);
        
        calibratedSpectra(t).meanAirTemperature=nanmean(vertcat(meteoData(rowInd).air_temperature))+273.15;
        calibratedSpectra(t).meanRelHumidity=nanmean(vertcat(meteoData(rowInd).rel_humidity));
        calibratedSpectra(t).meanAirPressure=nanmean(vertcat(meteoData(rowInd).air_pressure));
        calibratedSpectra(t).rainAccumulation=nansum(vertcat(meteoData(rowInd).rain_accumulation));
    end
else
    disp('format of meteo data changed before the 10th of August 2017') 
    
end
