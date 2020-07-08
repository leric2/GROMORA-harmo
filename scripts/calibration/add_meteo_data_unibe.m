function calibratedSpectra = add_meteo_data_unibe(calibrationTool,calibratedSpectra)
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


if datetime(str2num(calibrationTool.dateStr(1:4)),str2num(calibrationTool.dateStr(6:7)),str2num(calibrationTool.dateStr(9:10))) > datetime(2017,08,10)
    % First reading the Meteo dataset for this day
    dateStringMeteo=[calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)];
    meteoDataFile=[calibrationTool.meteoFolder 'meteo_' dateStringMeteo '.csv'];
    
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
        calibratedSpectra(t).meanRelHumidity=0.01*nanmean(vertcat(meteoData(rowInd).rel_humidity));
        calibratedSpectra(t).meanAirPressure=nanmean(vertcat(meteoData(rowInd).air_pressure));
        calibratedSpectra(t).rainAccumulation=nansum(vertcat(meteoData(rowInd).rain_accumulation));
    end
    
else
    disp('format of meteo data changed before the 10th of August 2017') 
    
end
