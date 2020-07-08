function meteoData = read_meteo_data_unibe(calibrationTool)
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
    
    meteoData(1).precipitation = meteoData(1).rain_accumulation;
    for i = 1:length(meteoData)
        meteoData(i).dateTime=datenum(meteoData(i).time,'yyyy-mm-ddTHH:MM:SS.FFFZ')-datenum(1970,1,1);
        
        % TODO Check units
        if i>1
            meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
        end
    end
else
    disp('format of meteo data changed before the 10th of August 2017') 
    
end
