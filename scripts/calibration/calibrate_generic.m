function [drift,calibratedSpectra] = calibrate_generic(rawSpectra,logFile,calibrationTool,calType)
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
%               |               - TCold
%               |               - Day, Month, Year
%               |               - timeZone
%               |               - elevationAngleHot, elevationAngleCold
%               |               - elevationAngleHotTol
%               |               - elevationAngleColdTol
%               |               - threshNumRawSpectraHot
%               |               - threshNumRawSpectraCold
%               |               - elevationAngleAntenna
%               |               - elevationAngleTolerance
%               |               - threshNumRawSpectraAnt
%               |               - numberOfChannels
%               |               - outlierDectectionType
%               |              
%               |         4. calType: Calibration type to do
%               |               1. standard: mean antenna vs mean hot/cold 
%               |                  spectra
%               |               2. debug: perform the standard one and:
%               |                  - Mean Up/Down: mean AT up, mean AT Down
%               |                  vs mean hot/cold
%               |                  - all cycles: all individual
%               |                  hot-cold-sky cycle BUT uncleaned for 
%               |                  potential spurious spectra ! Caution
%               |                  when comparing it to 'standard' mode if
%               |                  outliers in hot-cold are identified.
%               |                  
%               |
%               | OUTPUTS: 1. drift: structure containing mean values of all
%               |               channels for this day
%               |          2. calibratedSpectra: "structure array" of
%               |               calibrated data
%               |
% CALLS         | find_up_down_cycle()
%               |
%==========================================================================

% Calibration version
calibVersion = calibrationTool.calibrationVersion;

% CalibrationTime in Minute
calibTime=calibrationTool.calibrationTime;

% Extract all indices corresponding to hot, cold and sky observations in
% the log file.
if strcmp(calibrationTool.instrumentName,'mopi5')
    initialIndices={
        find(logFile.Position==calibrationTool.indiceHot & logFile.Measurement_NoiseDiode'==0);     % Hot
        find(logFile.Position==calibrationTool.indiceAntenna & logFile.Measurement_NoiseDiode'==0); % Antenna
        find(logFile.Position==calibrationTool.indiceCold & logFile.Measurement_NoiseDiode'==0);    % Cold
        };
else
    initialIndices={
        find(logFile.Position==calibrationTool.indiceHot & logFile.Tipping_Curve_active==0);      % Hot
        find(logFile.Position==calibrationTool.indiceAntenna & logFile.Tipping_Curve_active==0);  % Antenna
        find(logFile.Position==calibrationTool.indiceCold & logFile.Tipping_Curve_active==0);     % Cold
        };
end
    
if strcmp(calType,'debug')
    % Now we read complete half cycle Up (c-a-h) and Down (h-a-c) separately
    [firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(logFile,calibrationTool);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building the drift structure, inspired from Axel
% Variations of mean amplitude and TNoise with time
% need equal number of cycles !

for i=1:length(initialIndices)
    initialIndices{i}(initialIndices{i}>size(rawSpectra,1))=[];
    cn(i) = length(initialIndices{i}); 
end

% TODO in case the order of the cycle changes ! (for instance, if we have 2
% times more sky observations than cold or hot spectra, the following
% removal of indices fails !
driftIndices = initialIndices;
standardOrder = cn(1)/cn(2)>0.9 && cn(3)/cn(2)>0.9;
if ~standardOrder
    warning('order of the cycle is not standard, no drift quantities !')
else
    % we remove extra indices to keep only the smallest number of indices
    %between hot, cold and sky
    for i=1:length(initialIndices); 	driftIndices{i}=initialIndices{i}(1:min(cn)); end; 
end
dailyMeanTHot=mean(logFile.T_Hot_Absorber);
drift=struct();
% Drift always defined with hot indices
drift.t  = logFile.t(driftIndices{1});
drift.allDateTime = logFile.dateTime;
drift.cycleTime = diff(logFile.dateTime);
drift.dateTime  = logFile.dateTime(driftIndices{1});
drift.T  = logFile.T_Hot_Absorber(driftIndices{1});
if standardOrder
    drift.a(1,:) = mean(rawSpectra(driftIndices{1},:),2,'omitnan');
    drift.a(2,:) = mean(rawSpectra(driftIndices{2},:),2,'omitnan');
    drift.a(3,:) = mean(rawSpectra(driftIndices{3},:),2,'omitnan');
    drift.Y    = drift.a(1,:) ./  drift.a(3,:);
    drift.Tn   = (dailyMeanTHot - drift.Y*calibrationTool.TCold)./ (drift.Y-1);
    drift.Ta   = (drift.a(2,:) - drift.a(3,:)) ./ (drift.a(1,:) - drift.a(3,:)) *(dailyMeanTHot-calibrationTool.TCold) + calibrationTool.TCold;
else
    drift.a(1,:) = nan*ones(size(driftIndices{1}));
    drift.a(2,:) = nan*ones(size(driftIndices{1}));
    drift.a(3,:) = nan*ones(size(driftIndices{1}));
    drift.Y    = nan*ones(size(driftIndices{1}));
    drift.Tn   = nan*ones(size(driftIndices{1}));
    drift.Ta   = nan*ones(size(driftIndices{1}));
end
drift.TNoiseLog=logFile.FE_T_Sys(driftIndices{1});
drift.dailyMedianHotSpectra=median(rawSpectra(driftIndices{1},:),1,'omitnan');
drift.dailyStdHotSpectra=nanstd(rawSpectra(driftIndices{1},:));
drift.dailyMedianColdSpectra=median(rawSpectra(driftIndices{3},:),1,'omitnan');
drift.dailyStdColdSpectra=nanstd(rawSpectra(driftIndices{3},:));
drift.dailyMedianAntSpectra=median(rawSpectra(driftIndices{2},:),1,'omitnan');
drift.dailyStdAntSpectra=nanstd(rawSpectra(driftIndices{2},:));

drift.outlierCold = [];
drift.outlierHot = [];
drift.outlierSky = [];
drift.outlierDrift = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time interval for the calibration, splitting the indices for this day
% into the different calibration cycles

% Threshold for the separation, calibTime has to be in [min]
dt = hours(calibTime/60);

timeThresh = datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day, 'TimeZone', calibrationTool.timeZone):dt:datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day+1, 'TimeZone', calibrationTool.timeZone);

% Starting time for the complete cycle (will be the basis for separating
% the cycles according to calibrationTime). 
startingTimesHot=logFile.dateTime(initialIndices{1});
startingTimesAntennaAll=logFile.dateTime(initialIndices{2});
startingTimesCold=logFile.dateTime(initialIndices{3});

if strcmp(calType,'debug')
    startingTimesUp=logFile.dateTime(firstIndHalfUp);
    startingTimesDown=logFile.dateTime(firstIndHalfDown);
end
% Storing the indices specific to each calibration cycle in a new
% structure because by separating by time, we do not have the same
% number of individual cycle per calibration cycle
indices=struct();
for i = 1:length(timeThresh)-1
    % setting conditions based on the starting time of the individual cycles
    condAntenna=startingTimesAntennaAll>timeThresh(i) & startingTimesAntennaAll<timeThresh(i)+dt;   
    condHot=startingTimesHot>timeThresh(i) & startingTimesHot<timeThresh(i)+dt;
    condCold=startingTimesCold>timeThresh(i) & startingTimesCold<timeThresh(i)+dt;
    
    % saving the indices
    indices(i).validAntenna=initialIndices{2}(condAntenna);
    indices(i).validHot=initialIndices{1}(condHot);
    indices(i).validCold=initialIndices{3}(condCold);
    if strcmp(calType,'debug')
        % The same for Up and Down calibration
        condUp=startingTimesUp>timeThresh(i) & startingTimesUp<timeThresh(i)+dt;
        condDown=startingTimesDown>timeThresh(i) & startingTimesDown<timeThresh(i)+dt;
    
        indices(i).validColdStartUp=[firstIndHalfUp(condUp) firstIndHalfUp(condUp)+1 firstIndHalfUp(condUp)+2];
        indices(i).validHotStartDown=[firstIndHalfDown(condDown) firstIndHalfDown(condDown)+1 firstIndHalfDown(condDown)+2];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of the actual calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure containing all the information about the calibrated spectra
calibratedSpectra=struct();

% Number of calibration cycles for this day
nCalibrationCycles=length(indices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hot and Cold averaged spectra
for i=1:nCalibrationCycles

    calibratedSpectra(i).theoreticalStartTime=timeThresh(i);
    calibratedSpectra(i).theoreticalStartTimeStr=datestr(timeThresh(i),'yyyy_mm_dd_HH:MM:SS');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency vector
    if calibrationTool.IQProcessing
        calibratedSpectra(i).if = calibrationTool.samplingRateFFTS/2 * [-1:2/calibrationTool.numberOfChannels:1-2/calibrationTool.numberOfChannels];
        
        calibratedSpectra(i).observationFreq=calibrationTool.observationFreq;
        
        calibratedSpectra(i).LOFreqTot=calibrationTool.LOFreqTot;
        
        calibratedSpectra(i).freq=calibratedSpectra(i).if*1e6+calibratedSpectra(i).LOFreqTot;
        
        calibratedSpectra(i).df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    else
        calibratedSpectra(i).if = calibrationTool.samplingRateFFTS/2 * [0:1/calibrationTool.numberOfChannels:1-1/calibrationTool.numberOfChannels];
        
        calibratedSpectra(i).observationFreq=calibrationTool.observationFreq;
        
        calibratedSpectra(i).LOFreqTot=calibrationTool.LOFreqTot;
        
        calibratedSpectra(i).freq=calibratedSpectra(i).if*1e6+calibratedSpectra(i).LOFreqTot;
        
        calibratedSpectra(i).df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    end
    
    % hot and cold indices for this calib cycle
    ih=reshape(indices(i).validHot,1,[]);
    ic=reshape(indices(i).validCold,1,[]);

    initSizeHot=length(ih);
    initSizeCold=length(ic);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outlier detection for hot and cold spectra
    % Checking and removing any spurious angle for hot and cold   
    hotAngleOutlier=reshape((abs(logFile.Elevation_Angle(ih)-calibrationTool.elevationAngleHot) > calibrationTool.elevationAngleHotTol),[],1);
    coldAngleOutlier=reshape((abs(logFile.Elevation_Angle(ic)-calibrationTool.elevationAngleCold) > calibrationTool.elevationAngleColdTol),[],1);
       
    % identifying FFT spectrometer overloads 
    FFT_adc_overload_hot = reshape((logFile.FFT_adc_overload(ih) > calibrationTool.adcOverloadThresh),[],1);
    FFT_adc_overload_cold = reshape((logFile.FFT_adc_overload(ic) > calibrationTool.adcOverloadThresh),[],1);  
    
    % We compute the median of the calibration spectra as well as its std
    % deviation. We then remove the spectra which contains too many
    % channels that are beyond the (median +- n*stdDev)
    
    % Using the daily median and stddev to check for spurious hot and cold
    % spectra during this cycle:   
    medStdDevThreshHot=abs((rawSpectra(ih,:)-drift.dailyMedianHotSpectra))>calibrationTool.hotSpectraNumberOfStdDev*drift.dailyStdHotSpectra;
    medStdDevThreshCold=abs((rawSpectra(ic,:)-drift.dailyMedianColdSpectra))>calibrationTool.coldSpectraNumberOfStdDev*drift.dailyStdColdSpectra;
    
    outlierDetectHot = reshape(sum(medStdDevThreshHot,2)>calibrationTool.threshNumRawSpectraHot,[],1);
    outlierDetectCold = reshape(sum(medStdDevThreshCold,2)>calibrationTool.threshNumRawSpectraCold,[],1);
    
    % on a short scale, taking the median of the calibration cycle instead 
    % of daily (not used effectively):    
    medStdDevThreshHotShort = abs((rawSpectra(ih,:)-nanmedian(rawSpectra(ih,:))))>calibrationTool.hotSpectraNumberOfStdDev*nanstd(rawSpectra(ih,:));
    medStdDevThreshColdShort = abs((rawSpectra(ic,:)-nanmedian(rawSpectra(ic,:))))>calibrationTool.hotSpectraNumberOfStdDev*nanstd(rawSpectra(ic,:));    
    
    outlierDetectHotShort = reshape(sum(medStdDevThreshHotShort,2)>calibrationTool.threshNumRawSpectraHot,[],1);
    outlierDetectColdShort = reshape(sum(medStdDevThreshColdShort,2)>calibrationTool.threshNumRawSpectraCold,[],1);
    
    medianColdSpectra = nanmedian(rawSpectra(ic,:));
    medianHotSpectra = nanmedian(rawSpectra(ih,:));
    
    hotSpectra = rawSpectra(ih,:);
    coldSpectra = rawSpectra(ic,:);
    switch calibrationTool.outlierDectectionType
        case {'RFI_removal','RFI'}
                % RFI filtering (in progress)
    % problematic channel
    RFI_bad_channels = (calibratedSpectra(i).if>301) & (calibratedSpectra(i).if<499);
    % Using the cycle median and stddev to check for spurious sky
    % spectra during this cycle:


    
    mm = 100;
    %calibratedSpectra(i).ColdSpectraRFI = rawSpectra(ic,:)-drift.dailyMedianColdSpectra;
    %calibratedSpectra(i).ColdSpectraRFI = movmean(rawSpectra(ic,:)-drift.dailyMedianColdSpectra,mm,2);
    calibratedSpectra(i).ColdSpectraRFI = movmean(rawSpectra(ic,:)-medianColdSpectra,mm,2);
    calibratedSpectra(i).HotSpectraRFI = movmean(rawSpectra(ih,:)-medianHotSpectra,mm,2);
    
    %calibratedSpectra(i).ColdSpectraRFI2 =  calibratedSpectra(i).ColdSpectraRFI - median(calibratedSpectra(i).ColdSpectraRFI(RFI_bad_channels),2);
    %calibratedSpectra(i).ColdSpectraRFI2 =  calibratedSpectra(i).ColdSpectraRFI - movmean(calibratedSpectra(i).ColdSpectraRFI,10*mm,2);
    %diff_daily = movmean(drift.dailyMedianColdSpectra-median(drift.dailyMedianColdSpectra),mm,2);
    
    mean_smoothed_diff_to_median_cold = abs(median(calibratedSpectra(i).ColdSpectraRFI(:,RFI_bad_channels),2));
    mean_smoothed_diff_to_median_hot = abs(median(calibratedSpectra(i).HotSpectraRFI(:,RFI_bad_channels),2));

    mean_smoothed_diff_to_median_cold(mean_smoothed_diff_to_median_cold<0.5) = 1 ;
    mean_smoothed_diff_to_median_hot(mean_smoothed_diff_to_median_hot<0.5) = 1 ;

    calibratedSpectra(i).freqRFI = calibratedSpectra(i).if + 3300;
    outlierRFICold = reshape(sum(abs(calibratedSpectra(i).ColdSpectraRFI(:,RFI_bad_channels)) > 10*mean_smoothed_diff_to_median_cold,2) > 10,[],1);
    outlierRFIHot = reshape(sum(abs(calibratedSpectra(i).HotSpectraRFI(:,RFI_bad_channels)) > 10*mean_smoothed_diff_to_median_hot,2) > 10,[],1);
        otherwise
    end
    %outlierRFICold2 = reshape(sum(abs(ColdRFI_removed) > 10*mean_smoothed_diff_to_median_cold,2) > 10,[],1);
    % Depending on the defined outlier detection technique for this
    % calibration, we use the identified outliers to remove the spurious
    % hot and cold individual spectra 
    switch calibrationTool.outlierDectectionType
        case 'standard' 
            outlierHot = (outlierDetectHot | hotAngleOutlier | FFT_adc_overload_hot);
            outlierCold = (outlierDetectCold | coldAngleOutlier | FFT_adc_overload_cold);
        case 'RFI_removal'
            outlierHot = (outlierRFIHot | outlierDetectHot | hotAngleOutlier | FFT_adc_overload_hot);
            outlierCold = (outlierRFICold | outlierDetectCold | coldAngleOutlier | FFT_adc_overload_cold);
        case 'RFI'

            outlierHot = (outlierDetectHot | hotAngleOutlier | FFT_adc_overload_hot);
            outlierCold = (outlierDetectCold | coldAngleOutlier | FFT_adc_overload_cold);
            coldSpectra(find(abs(calibratedSpectra(i).ColdSpectraRFI) > 1*mean_smoothed_diff_to_median_cold))=nan;
            hotSpectra(find(abs(calibratedSpectra(i).HotSpectraRFI) > 1*mean_smoothed_diff_to_median_hot))=nan;

            rawSpectra(ih,:) = hotSpectra;
            rawSpectra(ic,:) = coldSpectra;

                calibratedSpectra(i).outlierDetectHotRFI = sum(outlierRFIHot);
    calibratedSpectra(i).outlierDetectColdRFI = sum(outlierRFICold);
            %RFI_cold_channels_2_remove = find(abs(movmean(rawSpectra(ic,RFI_bad_channels)-medianColdSpectra(RFI_bad_channels),100,2)) > 10*mean_smoothed_diff_to_median_cold);
            %RFI_hot_channels_2_remove=  find(abs(movmean(rawSpectra(ih,RFI_bad_channels)-medianHotSpectra(RFI_bad_channels),100,2)) > 10*mean_smoothed_diff_to_median_hot,);
            %[row, colum] = rawSpectra(abs(movmean(rawSpectra(ic,:)-medianColdSpectra,100,2)) > 10);
%             rawSpectra(abs(movmean(rawSpectra(ih,:)-medianHotSpectra,100,2)) > 10)=nan;
            %rawSpectra(ic,RFI_cold_channels_2_remove') = nan;
            %rawSpectra(ih,RFI_hot_channels_2_remove') = nan;
        case 'noFFT'
            outlierHot = (outlierDetectHot | hotAngleOutlier);
            outlierCold = (outlierDetectCold | coldAngleOutlier);
        case 'none'
            outlierHot = zeros(size(ih));
            outlierCold = zeros(size(ic));
    end
    
    calibratedSpectra(i).outlierDetectHot = sum(outlierDetectHot);
    calibratedSpectra(i).outlierDetectCold = sum(outlierDetectCold);
    calibratedSpectra(i).outlierDetectHotShort = sum(outlierDetectHotShort);
    calibratedSpectra(i).outlierDetectColdShort = sum(outlierDetectColdShort);
    calibratedSpectra(i).hotAngleOutlier = sum(hotAngleOutlier);
    calibratedSpectra(i).coldAngleOutlier = sum(coldAngleOutlier);
    calibratedSpectra(i).FFT_adc_overload_hot = sum(FFT_adc_overload_hot);
    calibratedSpectra(i).FFT_adc_overload_cold =sum(FFT_adc_overload_cold);
    


    % And save the outliers dates for plotting them later
    if sum(outlierHot)>0
        drift.outlierHot = [drift.outlierHot; reshape(logFile.dateTime(ih(outlierHot)),[],1)];
    end
    
    if sum(outlierCold)>0
        drift.outlierCold = [drift.outlierCold; reshape(logFile.dateTime(ic(outlierCold)),[],1)];
    end
    
    ih=ih(~outlierHot);
    ic=ic(~outlierCold);
    
    calibratedSpectra(i).spuriousHotSpectra=initSizeHot-length(ih);
    calibratedSpectra(i).spuriousColdSpectra=initSizeCold-length(ic); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % additionnal check with the drift structure (computing stdTNoise for
    % flagging later).
    % NOT USED NOW !
    % 
    % Here we only use the cleaned hot indices (that's the ones used in
    % drift datetime) to extract the drift noise receiver temperature,
    % its mean and its std dev during this calibration cycle.
    % TODO: simplify this if...end and remove unused stuff
    if ~isempty(ih)
        % Use drift structure for additionnal quality check
%         k=zeros(1,length(ih));
%         for a = 1:length(ih)
%             k(a) = find(drift.dateTime == logFile.dateTime(ih(a)));
%         end
    
        %Tn_drift_i=drift.Tn(k);
        %Ta_drift_i=drift.Ta(k);
    
        %outlierDrift = (abs(Tn_drift_i-median(Tn_drift_i))>3*std(Tn_drift_i) | abs(Ta_drift_i-median(Ta_drift_i))>4*std(Ta_drift_i))';
        
%         if sum(outlierDrift)>0
%             k(outlierDrift) = [];
%             drift.outlierDrift = [drift.outlierDrift; reshape(drift.dateTime(k(outlierDrift)),[],1)];
%             indHot=ih(outlierDrift);
%             indCold = ih(outlierDrift)-3;
%             for out = 1:length(indHot)
%                 ih(ih == indHot(out)) = [];
%                 ic(ic == indCold(out)) = [];
%             end
%         end

        % Find corresponding time stamps in the drift structure for this
        % cycle. 
        k = find((drift.dateTime > logFile.dateTime(ih(1)) & drift.dateTime < logFile.dateTime(ih(end))));
        if ~isempty(k)
        % Also used for stddev TNoise
        calibratedSpectra(i).TNoiseDrift=drift.Tn(k);
        calibratedSpectra(i).meanTNoiseDrift=nanmean(drift.Tn(k));
        calibratedSpectra(i).stdTNoise=nanstd(drift.Tn(k));
        else
           	calibratedSpectra(i).TNoiseDrift=NaN;
            calibratedSpectra(i).meanTNoiseDrift=NaN;
            calibratedSpectra(i).stdTNoise=NaN;
        end
    else
        calibratedSpectra(i).TNoiseDrift=NaN;
        calibratedSpectra(i).meanTNoiseDrift=NaN;
        calibratedSpectra(i).stdTNoise=NaN;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computing additional data from hot and cold spectra
    
    % Saving clean hot and cold indices for this cycle
    calibratedSpectra(i).hotInd=ih;
    calibratedSpectra(i).coldInd=ic;
    
    % Final mean hot and cold raw counts for this cycle:
    calibratedSpectra(i).meanHotSpectra=nanmean(rawSpectra(ih,:),1);
    calibratedSpectra(i).meanColdSpectra=nanmean(rawSpectra(ic,:),1);
    
    % Final std dev hot and cold raw counts for this cycle:
    calibratedSpectra(i).stdHotSpectra=nanstd(rawSpectra(ih,:),1);
    calibratedSpectra(i).stdColdSpectra=nanstd(rawSpectra(ic,:),1);
    
    % Hot temperature recorded during hot and cold spectra observations
    % (but always on the hot load of course)
    calibratedSpectra(i).THot=nanmean(logFile.T_Hot_Absorber([ih,ic]));
    calibratedSpectra(i).stdTHot=nanstd(logFile.T_Hot_Absorber([ih,ic]));
    
    % Computation of spectral Y-factor and noise receiver temperature for this cycle
    calibratedSpectra(i).Yspectral=calibratedSpectra(i).meanHotSpectra./calibratedSpectra(i).meanColdSpectra;
    
    % Noise Temperature Spectra (uncleaned for outliers in Y)
    calibratedSpectra(i).TN=(calibratedSpectra(i).THot - calibratedSpectra(i).Yspectral*calibrationTool.TCold)./ (calibratedSpectra(i).Yspectral -1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration of the sky (antenna) measurement with hot-cold calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select calibration type ('standard' or 'debug')
switch calType
    case 'standard'
        % Averaging hot and cold FFTS counts on the time interval. 
        for i=1:nCalibrationCycles
            calibratedSpectra(i).calibrationType=calType;
            calibratedSpectra(i).calibrationTime=calibTime;
            calibratedSpectra(i).calibrationVersion=calibVersion;
            
            % All antenna measurements during this cycle
            ia=reshape(indices(i).validAntenna,1,[]);
            
            % Antenna measurements inside a half cycle
            %iaUp=reshape(indices(i).validColdStartUp(2,:),1,[]);
            %iaDown=reshape(indices(i).validHotStartDown(2,:),1,[]);
            
            % Checking for NaN in the antenna spectra and keeping only complete
            % spectra for the calibration (rare ?):
            ia=ia(sum(isnan(rawSpectra(ia,:)),2)<1);
            %iaUp=iaUp(sum(isnan(rawSpectra(iaUp,:)),2)<1);
            %iaDown=iaDown(sum(isnan(rawSpectra(iaDown,:)),2)<1);
            
            % Saving all the indices for the Antenna
            calibratedSpectra(i).antennaInd=ia;
            %calibratedSpectra(i).antennaIndUp=iaUp;
            %calibratedSpectra(i).antennaIndDown=iaDown;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Outlier detection for the antenna individual spectra
            % Pointing angles:
            skyAngleCheck=abs(logFile.Elevation_Angle(calibratedSpectra(i).antennaInd)-calibrationTool.elevationAngleAntenna)>calibrationTool.elevationAngleTolerance;
            skyAngleCheck=reshape(skyAngleCheck,[],1);
            
            % FFT overload check:
            FFT_adc_overload_sky = reshape((logFile.FFT_adc_overload(ia) > calibrationTool.adcOverloadThresh),[],1);
            
            % Using the cycle median and stddev to check for spurious sky
            % spectra during this cycle:   
            medianSpectra = nanmedian(rawSpectra(ia,:));
            stdAntSpectra = nanstd(rawSpectra(ia,:));
            
            medStdDevThreshSky=abs((rawSpectra(ia,:)-medianSpectra))>calibrationTool.skySpectraNumberOfStdDev*stdAntSpectra;
            outlierDetectSky = reshape(sum(medStdDevThreshSky,2)>calibrationTool.threshNumRawSpectraAnt,[],1);
           

            % Depending on the defined outlier detection technique for this
            % calibration, we use the identified outliers to remove the spurious
            % sky individual spectra
            switch calibrationTool.outlierDectectionType
                case 'standard'
                    outlierSky = (outlierDetectSky | skyAngleCheck | FFT_adc_overload_sky);
                case 'RFI_removal'
                    % RFI filtering (in progress)
                    skySpectra = rawSpectra(ia,:);
                    % problematic channel
                    RFI_bad_channels = (calibratedSpectra(i).if>301) & (calibratedSpectra(i).if<499);
                    % Using the cycle median and stddev to check for spurious sky
                    % spectra during this cycle:

                    calibratedSpectra(i).AntSpectraRFI = movmean(rawSpectra(ia,:)-medianSpectra,100,2);
                    mean_smoothed_diff_to_median = abs(median(calibratedSpectra(i).AntSpectraRFI(:,RFI_bad_channels),2));
                    mean_smoothed_diff_to_median(mean_smoothed_diff_to_median<0.5) = 1 ;

                    %calibratedSpectra(i).freqRFI = calibratedSpectra(i).if(RFI_bad_channels) + 3300;
                    outlierRFI = reshape(sum(abs(calibratedSpectra(i).AntSpectraRFI(:,RFI_bad_channels)) > 10*mean_smoothed_diff_to_median,2) > 10,[],1);
                    outlierSky = (outlierRFI | outlierDetectSky | skyAngleCheck | FFT_adc_overload_sky);
                    calibratedSpectra(i).outlierRFISky=sum(outlierRFI);
                case 'RFI'
                    % RFI filtering (in progress)
                    skySpectra = rawSpectra(ia,:);
                    % problematic channel
                    RFI_bad_channels = (calibratedSpectra(i).if>301) & (calibratedSpectra(i).if<499);
                    % Using the cycle median and stddev to check for spurious sky
                    % spectra during this cycle:

                    calibratedSpectra(i).AntSpectraRFI = movmean(rawSpectra(ia,:)-medianSpectra,100,2);
                    mean_smoothed_diff_to_median = abs(median(calibratedSpectra(i).AntSpectraRFI(:,RFI_bad_channels),2));
                    mean_smoothed_diff_to_median(mean_smoothed_diff_to_median<0.5) = 1 ;

                    %calibratedSpectra(i).freqRFI = calibratedSpectra(i).if(RFI_bad_channels) + 3300;
                    outlierRFI = reshape(sum(abs(calibratedSpectra(i).AntSpectraRFI(:,RFI_bad_channels)) > 10*mean_smoothed_diff_to_median,2) > 10,[],1);
                    outlierSky = (outlierDetectSky | skyAngleCheck | FFT_adc_overload_sky);
                    skySpectra(find(abs(calibratedSpectra(i).AntSpectraRFI) > 10*mean_smoothed_diff_to_median))=nan;

                    rawSpectra(ia,:) = skySpectra;
                    calibratedSpectra(i).outlierRFISky=sum(outlierRFI);
                case 'noFFT'
                    outlierSky = (outlierDetectSky | skyAngleCheck);
                case 'none'
                    outlierSky = zeros(size(ia));
            end
            
            calibratedSpectra(i).outlierDetectSky=sum(outlierDetectSky);
            calibratedSpectra(i).skyAngleCheck=sum(skyAngleCheck);
            calibratedSpectra(i).FFT_adc_overload_sky=sum(FFT_adc_overload_sky);
            
            calibratedSpectra(i).antennaIndCleanAngle=ia(~outlierSky);
            calibratedSpectra(i).spuriousSkySpectra=sum(outlierSky);
            
            % save sky spurious spectra for later plotting
            if sum(outlierSky)>0
                drift.outlierSky = [drift.outlierSky; reshape(logFile.dateTime(ia(outlierSky)),[],1)];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Applying the calibration formula for this calibration cycle:

            % Clean sky FFT counts on this cycle:
            rsAntennaAll=rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:);
            
            % Mean clean sky FFT counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);
            calibratedSpectra(i).meanSkySpectra = rsAntenna;
            % Calibration
            
            % For all antenna measurement (to get the stddev of Tb on the
            % calibration cycle) against the mean of hot and cold spectra
            % for this calibration cycle.
                        
            % In terms of Intensity (same as RJE)
            IColdPlanck = planck_function(calibrationTool, calibrationTool.TCold, calibratedSpectra(i).freq);
            IHotPlanck = planck_function(calibrationTool, calibratedSpectra(i).THot, calibratedSpectra(i).freq);
            intensityAll = IColdPlanck + (IHotPlanck-IColdPlanck).*(rsAntennaAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);

            % In terms of RJE brightness temperature (same as Intensity)
            %TbColdRJ = planck_function(calibrationTool, calibrationTool.TCold, calibratedSpectra(i).freq)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*( calibratedSpectra(i).freq.^2));
            %TbHotRJ = planck_function(calibrationTool, calibratedSpectra(i).THot, calibratedSpectra(i).freq)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*( calibratedSpectra(i).freq.^2));
            
            if calibrationTool.savePlanckIntensity
                TbAll = (calibrationTool.h*calibratedSpectra(i).freq/calibrationTool.kb)./log((2*calibrationTool.h*calibratedSpectra(i).freq.^3)./(intensityAll.*calibrationTool.lightSpeed^2) + 1);
            else
                TbAll = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            end
            
            if size(std(TbAll),2) == calibrationTool.numberOfChannels
                calibratedSpectra(i).stdTb = std(TbAll);
            else
                calibratedSpectra(i).stdTb = -9999*ones(1,calibrationTool.numberOfChannels);
            end
            
            calibratedSpectra(i).TbPhysicalTemperature = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);

            % And with the mean sky spectrum of the cycle (it should be the
            % same as just taking the average of TbAll).
            if calibrationTool.savePlanckIntensity
                calibratedSpectra(i).intensityPlanck = IColdPlanck + (IHotPlanck-IColdPlanck).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
                calibratedSpectra(i).intensityPlanck(calibratedSpectra(i).intensityPlanck<=0) = NaN;
                calibratedSpectra(i).Tb = planck_Tb(calibrationTool, calibratedSpectra(i).intensityPlanck, calibratedSpectra(i).freq);
            else
                calibratedSpectra(i).Tb = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra); 
            end
            %calibratedSpectra(i).TbRJE = TbColdRJ + (TbHotRJ-TbColdRJ).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
        end
    case 'debug'
        % In this mode, we compute 3 types of calibration for debugging
        % purpose:
        % 1. "standard" one: mean antenna vs mean hot/cold spectra
        % 2. "Mean Up/Down": mean AT up, mean AT Down vs mean hot/cold
        % 3. "all cycles": all individual cycle vs mean hot/cold,
        % uncleaned for pointing angle problems.
        %
        %
        % Averaging hot and cold FFTS counts on the time interval.
        % Keeping each individual cycle for later analysis but computing as
        % well the avg calibrated spectra.
        
        for i=1:nCalibrationCycles
            calibratedSpectra(i).calibrationVersion=calibVersion;
            calibratedSpectra(i).calibrationTime=calibTime;
            
            % All antenna measurements
            ia=reshape(indices(i).validAntenna,1,[]);
                     
            % Sky measurements inside a half cycle
            iaUp=reshape(indices(i).validColdStartUp(:,2),1,[]);
            iaDown=reshape(indices(i).validHotStartDown(:,2),1,[]);
            
            % Checking for NaN in the antenna spectra and keeping only complete
            % spectra for the calibration:
            ia=ia(sum(isnan(rawSpectra(ia,:)),2)<1);
            iaUp=iaUp(sum(isnan(rawSpectra(iaUp,:)),2)<1);
            iaDown=iaDown(sum(isnan(rawSpectra(iaDown,:)),2)<1);
            
            % Saving the indices for the Antenna
            calibratedSpectra(i).antennaInd=ia;
            calibratedSpectra(i).antennaIndUp=iaUp;
            calibratedSpectra(i).antennaIndDown=iaDown;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Outlier detection for the antenna individual spectra
            % Pointing angles:
            skyAngleCheck=abs(logFile.Elevation_Angle(calibratedSpectra(i).antennaInd)-calibrationTool.elevationAngleAntenna)>calibrationTool.elevationAngleTolerance;
            skyAngleCheck=reshape(skyAngleCheck,[],1);

            % FFT overload check:
            FFT_adc_overload_sky = reshape(~(logFile.FFT_adc_overload(ia) == 0),[],1);
            
            % Using the cycle median and stddev to check for spurious sky
            % spectra during this cycle:   
            medianSpectra = nanmedian(rawSpectra(ia,:));
            stdAntSpectra = nanstd(rawSpectra(ia,:));
            
            medStdDevThreshSky=abs((rawSpectra(ia,:)-medianSpectra))>calibrationTool.skySpectraNumberOfStdDev*stdAntSpectra;
            outlierDetectSky = reshape(sum(medStdDevThreshSky,2)>calibrationTool.threshNumRawSpectraAnt,[],1);
            
            % Depending on the defined outlier detection technique for this
            % calibration, we use the identified outliers to remove the spurious
            % sky individual spectra
            switch calibrationTool.outlierDectectionType
                case 'standard'
                    outlierSky = (outlierDetectSky | skyAngleCheck | FFT_adc_overload_sky);
                case 'noFFT'
                    outlierSky = (outlierDetectSky | skyAngleCheck);
                case 'none'
                    outlierSky = zeros(size(ia));
            end
            
            calibratedSpectra(i).antennaIndCleanAngle=ia(~outlierSky);
            calibratedSpectra(i).spuriousSkySpectra=sum(outlierSky);
            
            if sum(outlierSky)>0
                drift.outlierSky = [drift.outlierSky; logFile.dateTime(ia(outlierSky))];
            end
            
            % to remove also in Up Down
            ind2remove = ia(outlierSky);
            for l=1:length(ind2remove)              
                iaUp(iaUp == ind2remove(l)) = [];
                iaDown(iaDown == ind2remove(l)) = [];
            end
            
            calibratedSpectra(i).antennaIndCleanUpAngle=iaUp;
            calibratedSpectra(i).antennaIndCleanDownAngle=iaDown;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType=calType;
            
            % All clean Antenna counts on this cycle to compute the std Tb
            rsAntennaAllClean=rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:);
            
            TbAllCleanAngle = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaAllClean-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            calibratedSpectra(i).stdTb=std(TbAllCleanAngle);
            
            % Mean Antenna counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);
            calibratedSpectra(i).meanSkySpectra = rsAntenna;
            % Calibration 1. "standard"
            calibratedSpectra(i).Tb = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra); 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration separately (Mean Up and Down) for this calibration cycle:
            % Mean Antenna Up/Down for this cycle
            rsAntennaUp=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanUpAngle,:),1);
            rsAntennaDown=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanDownAngle,:),1);
            
            % Calibration 2. "Mean Up/Down"
            calibratedSpectra(i).TbUp = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaUp-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            calibratedSpectra(i).TbDown = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaDown-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calibration of all individual cycles separately. 
            % Saving the cycle by cycle calibration (separated for clarity)
            
            % For now, we do not keep the difference between up/dow cycle.
            % If one exist, it should be visible from "Mean Up/Down" type
            % of calibration (assuming Up/Down do not influence the mean
            % hot and cold spectra
            
            % Collecting all individual cycles around the clean sky
            % spectra. !! We do not recheck for cold and hot spectra
            % quality at that point !!
            rsAntennaUp = rawSpectra(calibratedSpectra(i).antennaIndCleanUpAngle,:);
            rsColdUp = rawSpectra(calibratedSpectra(i).antennaIndCleanUpAngle-1,:);
            rsHotUp = rawSpectra(calibratedSpectra(i).antennaIndCleanUpAngle+1,:);
            
            rsAntennaDown=rawSpectra(calibratedSpectra(i).antennaIndCleanDownAngle,:);
            rsColdDown = rawSpectra(calibratedSpectra(i).antennaIndCleanDownAngle+1,:);
            rsHotDown = rawSpectra(calibratedSpectra(i).antennaIndCleanDownAngle-1,:);
            
            % Calibration 3. "all cycles"
            calibratedSpectra(i).TbAllUp = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaUp-rsColdUp)./(rsHotUp-rsColdUp);
            calibratedSpectra(i).TbAllDown = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaDown-rsColdDown)./(rsHotDown-rsColdDown);
            calibratedSpectra(i).TbAll = [calibratedSpectra(i).TbAllUp ; calibratedSpectra(i).TbAllDown];
            
            % TODO:            
            % To check the difference if we do the averaging after the
            % calibration of individual cycle:
            % Compute mean cold and hot spectra from ALL individual cycles
            % (also the uncleaned ones) and compare it to the avg TbAll 
            % from above --> Use the 'none' outliers removal techniques ?!
            calibratedSpectra(i).TbFromAllAvgAfter=nanmean(calibratedSpectra(i).TbAll);
            
            % Sky outliers
            calibratedSpectra(i).skyFlag=outlierSky;
        end
end
end

