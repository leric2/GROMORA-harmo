function calibratedSpectra = check_calibrated_miawara_LN2(standardLog,calibrationTool,calibratedSpectra)
%==========================================================================
% NAME      | check_calibrated_generic.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2020
%           |
% ABSTRACT  | Function performing some checks on the calibration and
%           | adding meta information to the calibratedSpectra structure.
%           | It builts also the flags vector for the level1a data for each
%           | calibration cycle and add every new information to the
%           | calibrated spectra structure (IN/OUT).
%           | 
%           |
% ARGUMENTS | INPUTS:   - standardLog: harmonized GROSOM log file 
%           |           - calibrationTool
%           |           - calibratedSpectra
%           |
%           |
%           | OUTPUTS: - calibratedSpectra
%           |
% CALLS     | 
%           | 
%==========================================================================
% Checking all calibration cycle

for i = 1:size(calibratedSpectra,2)
    f = load(calibrationTool.channel_freqs);
    calibratedSpectra(i).freq = f(:,2);
    

        f  = load(calibrationTool.channel_freqs);
        N  = calibrationTool.numberOfChannels;
        calibratedSpectra(i).observationFreq=f;
        calibratedSpectra(i).freq = f;%N/2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        calibratedSpectra(i).year=year(datenum(calibrationTool.dateStr));
        calibratedSpectra(i).month=month(datenum(calibrationTool.dateStr));
        calibratedSpectra(i).day=day(datenum(calibrationTool.dateStr));

        calibratedSpectra(i).timeMin = datestr(calibratedSpectra(i).theoreticalStartTime,'YYYY_mm_dd_HH:MM:SS');
        calibratedSpectra(i).timeMax = datestr(calibratedSpectra(i).theoreticalStartTime + minutes(calibratedSpectra(i).calibrationTime),'YYYY_mm_dd_HH:MM:SS');

        calibratedSpectra(i).timeMin=datenum(calibratedSpectra(i).timeMin,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
        calibratedSpectra(i).timeMax=datenum(calibratedSpectra(i).timeMax,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);


        if standardLog.Month(1) < 10
            m = ['0' num2str(standardLog.Month(1))];
        else
            m = num2str(standardLog.Month(1));
        end

        if standardLog.Day(1) < 10
            d = ['0' num2str(standardLog.Day(1))];
        else
            d = num2str(standardLog.Day(1));
        end

        calibratedSpectra(i).date=calibrationTool.dateStr;   %%[num2str(standardLog.Year(1)) '_' m '_' d];

       
        meanDatetime=[calibratedSpectra(i).date '_' datestr(mean(standardLog.t(:))/24,'HH:MM:SS')];
        calibratedSpectra(i).meanDatetime=datenum(meanDatetime,'YYYY_mm_dd_HH:MM:SS')-datenum(1970,1,1);
        calibratedSpectra(i).timeOfDay=mean(standardLog.t(:));

       
        idx = find(calibratedSpectra(i).freq(:,2)>=2.1786e+10 & calibratedSpectra(i).freq(:,2)<=2.2684e+10);

%         m            = length(idx);
%         [p,s,mu]     = polyfit(1:m,calibratedSpectra(i).TCold(idx),1);
%         level1_sigma = calibratedSpectra(i).TCold(idx)-polyval(p,1:m,[],mu);
%         calibratedSpectra(i).sigma        = std(level1_sigma);
        
        % Error vector description:
        calibratedSpectra(i).errorVectorDescription=[...
            'sufficientNumberOfIndices',...
            'systemTemperatureOK',...
            'LN2SensorsOK',...
            'LN2LevelOK',...
            'hotLoadOK',...
            'FFT_adc_overload_OK'];

      end
end

