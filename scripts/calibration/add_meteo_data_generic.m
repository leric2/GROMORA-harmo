function calibratedSpectra = add_meteo_data_generic(calibrationTool, meteoData, calibratedSpectra)
%==========================================================================
% NAME          | add_meteo_data_generic
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Function to add the meteo value to the calibrated spectra
%               | structure. It does a simple averaged on the values from
%               | meteoData depeding on the calibrationTime defined.
%               |
%               |
% ARGUMENTS     | INPUTS:  1. calibrationTool:
%               |            - referenceTime
%               |          2. meteoData: standard structure containing the
%               |            meteo data.
%               |          3. calibratedSpectra
%               |
%               | OUTPUTS: 1. calibratedSpectra with meteoData
%               |
% COMMENTS      | Called only from run_integration() before performing the 
%               | integration because meteoData are saved in full in
%               | level1a
%               |
%==========================================================================
%%%%%% Storing values into calibratedSpectra
for t=1:length(calibratedSpectra)
    %start=datetime(calibratedSpectra(t).timeMin+calibrationTool.referenceTime,'ConvertFrom','datenum');
    %stop=datenum(datetime(calibratedSpectra(t).timeMin+calibrationTool.referenceTime,'ConvertFrom','datenum')+seconds(calibratedSpectra(t).calibrationTime))-datenum(1970,1,1);
    timeMin = datetime(calibratedSpectra(t).time_min+calibrationTool.referenceTime,'ConvertFrom','datenum', 'TimeZone',calibrationTool.timeZone);
    stop = timeMin+seconds(calibratedSpectra(t).calibration_time);
    
    % Selecting the interesting values for each calibration cycle:
    rowInd=([meteoData.dateTime]>=timeMin & [meteoData.dateTime]<=stop);
    
    calibratedSpectra(t).mean_air_temperature=nanmean(vertcat(meteoData.air_temperature(rowInd)));
    calibratedSpectra(t).mean_relative_humidity=0.01*nanmean(vertcat(meteoData.relative_humidity(rowInd)));
    calibratedSpectra(t).mean_air_pressure=nanmean(vertcat(meteoData.air_pressure(rowInd)));
    calibratedSpectra(t).rain_accumulation=sum(vertcat(meteoData.precipitation(rowInd)),'omitnan');
end
end
