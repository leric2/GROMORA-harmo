function calibratedSpectra = add_meteo_data_generic(calibrationTool, meteoData, calibratedSpectra)
%==========================================================================
% NAME          | add_meteo_data_generic
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014 // modified 01.2020
%               |
% ABSTRACT      | Function to add the meteo value to the calibrated spectra
%               | 
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

%%%%%% Storing values into calibratedSpectra
for t=1:length(calibratedSpectra)
    %start=datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum');
    %stop=datenum(datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum')+seconds(calibratedSpectra(t).calibrationTime))-datenum(1970,1,1);
    timeMin = datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum');
    stop = timeMin+seconds(calibratedSpectra(t).calibrationTime);
    % Selecting the interesting values for each calibration cycle:
    rowInd=([meteoData.dateTime]>=timeMin & [meteoData.dateTime]<=stop);
    
    calibratedSpectra(t).meanAirTemperature=nanmean(vertcat(meteoData.air_temperature(rowInd)));
    calibratedSpectra(t).meanRelHumidity=0.01*nanmean(vertcat(meteoData.rel_humidity(rowInd)));
    calibratedSpectra(t).meanAirPressure=nanmean(vertcat(meteoData.air_pressure(rowInd)));
    calibratedSpectra(t).rainAccumulation=nansum(vertcat(meteoData.precipitation(rowInd)));
end
end
