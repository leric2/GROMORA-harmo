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
    timeMin = datetime(calibratedSpectra(t).time_min+datenum(1970,1,1),'ConvertFrom','datenum');
    stop = timeMin+seconds(calibratedSpectra(t).calibration_time);
    % Selecting the interesting values for each calibration cycle:
    rowInd=([meteoData.dateTime]>=timeMin & [meteoData.dateTime]<=stop);
    
    calibratedSpectra(t).mean_air_temperature=nanmean(vertcat(meteoData.air_temperature(rowInd)));
    calibratedSpectra(t).mean_relative_humidity=0.01*nanmean(vertcat(meteoData.relative_humidity(rowInd)));
    calibratedSpectra(t).mean_air_pressure=nanmean(vertcat(meteoData.air_pressure(rowInd)));
    calibratedSpectra(t).rain_accumulation=nansum(vertcat(meteoData.precipitation(rowInd)));
end
end
