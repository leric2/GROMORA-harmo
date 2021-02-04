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
tippingCurveSkyInd=find(logFile.Tipping_Curve_active & logFile.Position == calibrationTool.indiceTC); 

% the tipping curve is calculated for the average of a certain frequency
% range

% f = load(calibrationTool.channel_freqs);
N = calibrationTool.numberOfChannels;
% freq     = interp1(f(:,1)',f(:,2)',1:N/2);
% idx_freq = find (freq > 22.135e9 & freq < 22.335e9);

skipChannels = calibrationTool.TC.skipFraction * N;
        
lower = int16(skipChannels);
upper = int16(skipChannels) + calibrationTool.TC.numberOfChannelsTropCorr;

if strcmp(calibrationTool.TC.useWings,'both')
    idxFreqTC = [lower:upper,N-upper:N-lower];
elseif strcmp(calibrationTool.TC.useWings,'left')
    idxFreqTC = lower:upper;
elseif strcmp(calibrationTool.TC.useWings,'right')
    idxFreqTC = N-upper:N-lower;
else
    idxFreqTC = calibrationTool.TC.tippingCurveChannels;
end

%TC_data.tippingCurveMeanRawCounts=mean(rawSpectra(tippingCurveInd,idxFreqTC),2);
%TC_data.datetime = logFile.dateTime(tippingCurveInd);

%if length(tippingCurveInd) == length(tippingCurveSkyInd)
%    lastTipAngle = find(diff( tippingCurveSkyInd) > 10);
%    nTippingCurve = length(lastTipAngle);
%    if all(diff(lastTipAngle) == calibrationTool.tippingSize)
%        effTippingSize = calibrationTool.tippingSize;
%    else
%        effTippingSize = median(diff(lastTipAngle));
%    end
%    firstTipAngle=lastTipAngle-effTippingSize+1;
%else
    lastTipAngle = find(diff( tippingCurveInd) > 10);
    nTippingCurve = length(lastTipAngle);
    if all(diff(lastTipAngle) == calibrationTool.tippingSize)
       effTippingSize = calibrationTool.tippingSize;
    else
       effTippingSize = median(diff(lastTipAngle));
    end
    firstTipAngle=lastTipAngle-effTippingSize+1;
%end

for i = 1:nTippingCurve
    TC_data(i).position = logFile.Position(tippingCurveInd(firstTipAngle(i):lastTipAngle(i)))';
    TC_data(i).angle = logFile.Elevation_Angle(tippingCurveInd(firstTipAngle(i):lastTipAngle(i)))';
    TC_data(i).THot_all = logFile.T_Hot_Absorber(tippingCurveInd(firstTipAngle(i):lastTipAngle(i)))';
    TC_data(i).tippingCurveRawCounts = rawSpectra(tippingCurveInd(firstTipAngle(i):lastTipAngle(i)),idxFreqTC);
    TC_data(i).tippingCurveMeanRawCounts = nanmean(rawSpectra(tippingCurveInd(firstTipAngle(i):lastTipAngle(i)),idxFreqTC),2)';
    TC_data(i).datetime = logFile.dateTime(tippingCurveInd(firstTipAngle(i):lastTipAngle(i)))';
    TC_data(i).channels = idxFreqTC;
    TC_data(i).datetime.TimeZone='Z';
      
    if length(tippingCurveInd) ~= length(tippingCurveSkyInd)
        TC_data(i).cold = mean(TC_data(i).tippingCurveMeanRawCounts(TC_data(i).position == calibrationTool.indiceCold));
        TC_data(i).cold_spectraTC = TC_data(i).tippingCurveRawCounts(TC_data(i).position == calibrationTool.indiceCold,:);
        TC_data(i).hot = mean(TC_data(i).tippingCurveMeanRawCounts(TC_data(i).position == calibrationTool.indiceHot));
        TC_data(i).hot_spectraTC = TC_data(i).tippingCurveRawCounts(TC_data(i).position == calibrationTool.indiceHot,:);
        TC_data(i).THot = mean(TC_data(i).THot_all(TC_data(i).position == calibrationTool.indiceHot));
        TC_data(i).sky = TC_data(i).tippingCurveMeanRawCounts(TC_data(i).position == calibrationTool.indiceTC);
        TC_data(i).sky_spectra = TC_data(i).tippingCurveRawCounts(TC_data(i).position == calibrationTool.indiceTC,:);
        TC_data(i).skyAngle = TC_data(i).angle(TC_data(i).position == calibrationTool.indiceTC);
        TC_data(i).skyMeanDatetime = mean(TC_data(i).datetime(TC_data(i).position == calibrationTool.indiceTC));
        TC_data(i).meanDateNum = datenum(TC_data(i).skyMeanDatetime)-calibrationTool.referenceTime;
    else 
        TC_data(i).sky = TC_data(i).tippingCurveMeanRawCounts;
        TC_data(i).sky_spectra = TC_data(i).tippingCurveRawCounts;
        TC_data(i).skyAngle = TC_data(i).angle;
        TC_data(i).skyMeanDatetime = mean(TC_data(i).datetime);
        TC_data(i).meanDateNum = datenum(TC_data(i).skyMeanDatetime)-calibrationTool.referenceTime;
    end
end
end  
