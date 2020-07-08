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
    stop=datenum(datetime(calibratedSpectra(t).timeMin+datenum(1970,1,1),'ConvertFrom','datenum')+seconds(calibratedSpectra(t).calibrationTime))-datenum(1970,1,1);
    
    % Selecting the interesting values for each calibration cycle:
    rowInd=([meteoData.dateTime]>=calibratedSpectra(t).timeMin & [meteoData.dateTime]<=stop);
    
    calibratedSpectra(t).meanAirTemperature=nanmean(vertcat(meteoData(rowInd).air_temperature))+273.15;
    calibratedSpectra(t).meanRelHumidity=0.01*nanmean(vertcat(meteoData(rowInd).rel_humidity));
    calibratedSpectra(t).meanAirPressure=nanmean(vertcat(meteoData(rowInd).air_pressure));
    calibratedSpectra(t).rainAccumulation=nansum(vertcat(meteoData(rowInd).rain_accumulation));
end
end
