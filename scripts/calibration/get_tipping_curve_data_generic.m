function TC_data = get_tipping_curve_data_generic(rawSpectra, logFile, calibrationTool)
%==========================================================================
% NAME          | calibrate_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Reading and saving tipping curve calibration data within
%               | dedicated structure.
%               | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: 1. rawSpectra: Matrix of raw data (1 line per
%               |           cycle, channels as columns.
%               |         2. logFile: standardized log file
%               |         3. calibrationTool:
%               |               - numberOfChannels
%               |               - TC
%               |               - indiceHot, indiceAntenna, indiceCold,
%               |                 indiceTC
%               |               - referenceTime
%               |               - tippingSize
%               |              
%               |
%               | OUTPUTS: 1. TC_data: structure array containing all
%               |               tipping curve data for this day
%               |
%==========================================================================

% find TC cycles
tippingCurveInd=find(logFile.Tipping_Curve_active);      % all TC
tippingCurveSkyInd=find(logFile.Tipping_Curve_active & logFile.Position == calibrationTool.indiceTC); 

N = calibrationTool.numberOfChannels;

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

lastTipAngle = find(diff( tippingCurveInd) > 10);
nTippingCurve = length(lastTipAngle);
if all(diff(lastTipAngle) == calibrationTool.tippingSize)
    effTippingSize = calibrationTool.tippingSize;
else
    effTippingSize = median(diff(lastTipAngle));
end
firstTipAngle=lastTipAngle-effTippingSize+1;

% removing the first TC when it is not completed during this day
nTippingCurve = nTippingCurve - sum(firstTipAngle < 0);
lastTipAngle(firstTipAngle < 0) = [];
firstTipAngle(firstTipAngle < 0) = [];

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

    if i==1
        lenTC = length(TC_data(i).skyAngle);
    else
        if length(TC_data(i).skyAngle) ~= lenTC
            %TC_data(i) = [];
            % If the number of angles is not constant, we just don't do any
            % tipping curve (problem when saving the data). 
            error('Problem with tipping curve for this day, deactivating')
        end
    end
end
end  

