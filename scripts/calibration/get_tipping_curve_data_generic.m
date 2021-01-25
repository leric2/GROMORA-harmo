function TC_data = get_tipping_curve_data_generic(rawSpectra, logFile, calibrationTool)
%==========================================================================
% NAME          | calibrate_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Doing a hot-cold calibration for MW radiometer.
%               | This function has different modes to perform the
%               | calibration: 'standard' and 'debug' (see below) and
%               | different options to remove (or not) outliers
%               | (individual spurious spectrum) before averaging them
%               | together on the calibration time.
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: 1. rawSpectra: Matrix of raw data (1 line per
%               |           cycle, channels as columns.
%               |         2. logFile: standardized log file
%               |         3. calibrationTool:
%               |               - calibrationVersion
%               |               - calibrationTime
%               |               - instrumentName
%               |               - indiceHot, indiceAntenna, indiceCold

%               |              
%               |         4. calType: Calibration type to do
%               |               1. standard: mean antenna vs mean hot/cold 
%               |                  spectra
%               |               2. debug: perform the standard one and:
%               |
%               | OUTPUTS: 1. drift: structure containing mean values of all
%               |               channels for this day
%               |          2. calibratedSpectra: "structure array" of
%               |               calibrated data
%               |
% CALLS         | find_up_down_cycle()
%               |
%==========================================================================
% idx hot

%idx_hot = find(logFile.Mirror_pos == 1);

tippingCurveInd=find(logFile.Tipping_Curve_active);      % all TC

% the tipping curve is calculated for the average of a certain frequency
% range

% f = load(calibrationTool.channel_freqs);
% N        = calibrationTool.numberOfChannels;
% freq     = interp1(f(:,1)',f(:,2)',1:N/2);
% idx_freq = find (freq > 22.135e9 & freq < 22.335e9);

idxFreqTC = calibrationTool.tippingCurveChannels;

%TC_data.tippingCurveMeanRawCounts=mean(rawSpectra(tippingCurveInd,idxFreqTC),2);
%TC_data.datetime = logFile.dateTime(tippingCurveInd);

TippingCounter = 1;
for i = 1:length(tippingCurveInd)/calibrationTool.tippingSize
    TC_data(i).position = logFile.Position(tippingCurveInd(calibrationTool.tippingSize*(i-1)+1:calibrationTool.tippingSize*i));
    TC_data(i).angle = logFile.Elevation_Angle(tippingCurveInd(calibrationTool.tippingSize*(i-1)+1:calibrationTool.tippingSize*i));
    TC_data(i).THot = logFile.T_Hot_Absorber(tippingCurveInd(calibrationTool.tippingSize*(i-1)+1:calibrationTool.tippingSize*i));
    TC_data(i).tippingCurveMeanRawCounts = mean(rawSpectra(tippingCurveInd(calibrationTool.tippingSize*(i-1)+1:calibrationTool.tippingSize*i),idxFreqTC),2);
    TC_data(i).datetime = logFile.dateTime(tippingCurveInd(calibrationTool.tippingSize*(i-1)+1:calibrationTool.tippingSize*i));
    TC_data(i).datetime.TimeZone='Z';
      
    if strcmp(calibrationTool.TC_type, 'SkyLoads')
        TC_data(i).cold = mean(TC_data(i).tippingCurveMeanRawCounts(TC_data(i).position == calibrationTool.indiceCold));
        TC_data(i).hot = mean(TC_data(i).tippingCurveMeanRawCounts(TC_data(i).position == calibrationTool.indiceHot));
        TC_data(i).THot = mean(TC_data(i).THot(TC_data(i).position == calibrationTool.indiceHot));
        TC_data(i).sky = TC_data(i).tippingCurveMeanRawCounts(TC_data(i).position == calibrationTool.indiceTC)';
        TC_data(i).skyAngle = TC_data(i).angle(TC_data(i).position == calibrationTool.indiceTC)';
        TC_data(i).skyMeanDatetime = mean(TC_data(i).datetime(TC_data(i).position == calibrationTool.indiceTC));
        TC_data(i).meanDateNum = datenum(TC_data(i).skyMeanDatetime)-calibrationTool.referenceTime;
    elseif strcmp(calibrationTool.TC_type, 'onlySkyObs')
        TC_data(i).sky = TC_data(i).tippingCurveMeanRawCounts;
        TC_data(i).skyAngle = TC_data(i).angle;
        TC_data(i).skyMeanDatetime = mean(TC_data(i).datetime);
        TC_data(i).meanDateNum = datenum(TC_data(i).skyMeanDatetime)-calibrationTool.referenceTime;
    else
        error('TC type not recognized')
    end
end
end  

