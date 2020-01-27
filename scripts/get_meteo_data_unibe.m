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


if correctedSpectra.meanTime > (datenum(2017,08,10)-datenum(1970,1,1));
    dateStringMeteo=[retrievalTool.dateStr(1:4) '-' retrievalTool.dateStr(6:7) '-' retrievalTool.dateStr(9:10)];
    meteoDataFile=[retrievalTool.meteoFolder 'meteo_' dateStringMeteo '.csv'];
    
    % Reading meteo data and transforming it into matlab structure
    T=readtable(meteoDataFile);
    meteoData=table2struct(T);
    
    for i = 1:length(meteoData)
        meteoData(i).dateTime=datenum(meteoData(i).time,'yyyy-mm-ddTHH:MM:SS.FFFZ')-datenum(1970,1,1);
    end
    
    airTempVec=[meteoData.air_temperature];
    relHumVect=[meteoData.rel_humidity];
    
    meanTemp=ones(length(correctedSpectra.meanTime),1)*NaN;
    meanHumi=ones(length(correctedSpectra.meanTime),1)*NaN;
    
    for i = 1:length(correctedSpectra.meanTime)
        rowInd=([meteoData.dateTime]>correctedSpectra.firstSkyTime(i)) & ([meteoData.dateTime]<correctedSpectra.lastSkyTime(i));
        
        meanTemp(i)=nanmean(airTempVec(rowInd));
        meanHumi(i)=nanmean(relHumVect(rowInd));
    end
    
    correctedSpectra.meanTemperature=meanTemp;
    correctedSpectra.meanRelHumidity=meanHumi;
else
    disp('format of meteo data changed before the 10th of August 2017') 
end