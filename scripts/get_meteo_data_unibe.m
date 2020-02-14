function correctedSpectra = get_meteo_data_unibe(correctedSpectra,retrievalTool)
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

if datetime(correctedSpectra(1).year,correctedSpectra(1).month,correctedSpectra(1).day) > datetime(2017,08,10)
    % First reading the Meteo dataset for this day
    dateStringMeteo=[retrievalTool.dateStr(1:4) '-' retrievalTool.dateStr(6:7) '-' retrievalTool.dateStr(9:10)];
    meteoDataFile=[retrievalTool.meteoFolder 'meteo_' dateStringMeteo '.csv'];
    
    % Transforming it into matlab structure
    T=readtable(meteoDataFile);
    meteoData=table2struct(T);
    
    for i = 1:length(meteoData)
        meteoData(i).dateTime=datenum(meteoData(i).time,'yyyy-mm-ddTHH:MM:SS.FFFZ')-datenum(1970,1,1);
    end
    
    % Storing interesting variable into vectors
    airTempVec=[meteoData.air_temperature];
    relHumVect=[meteoData.rel_humidity];
    airPVec=[meteoData.air_pressure];
    
    for t=1:length(correctedSpectra)
        
        % Selecting the interesting values for each calibration cycle:
        rowInd=([meteoData.dateTime]>correctedSpectra(t).firstSkyTime) & ([meteoData.dateTime]<correctedSpectra(t).lastSkyTime);

        correctedSpectra(t).meanAirTemperature=nanmean(airTempVec(rowInd))+273.15;
        correctedSpectra(t).meanRelHumidity=nanmean(relHumVect(rowInd));
        correctedSpectra(t).meanAirPressure=nanmean(airPVec(rowInd));
    end
else
    disp('format of meteo data changed before the 10th of August 2017') 
    
end
