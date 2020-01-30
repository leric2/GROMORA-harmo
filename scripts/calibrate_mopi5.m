function calibratedSpectra = calibrate_mopi5(rawSpectra,log,retrievalTool,TCold,calType)
%==========================================================================
% NAME          | calibrate_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Doing a hot-cold calibration for MW radiometer.
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

% Calibration version
calibVersion='1.0.0';

% CalibrationTime in Minute
calibTime=retrievalTool.calibrationTime;

% Finding all indices of each type for this day that are not part of a
% tipping curve calibration
hotIndices=find(log.Position==retrievalTool.indiceHot & log.Measurement_NoiseDiode'==0);
antennaIndices=find(log.Position==retrievalTool.indiceAntenna & log.Measurement_NoiseDiode'==0);
coldIndices=find(log.Position==retrievalTool.indiceCold & log.Measurement_NoiseDiode'==0);

% % Indices for C-H-A with ND off/on
% C={
% 	find(log.Measurement_Target==2 &  log.Measurement_NoiseDiode==0); % Cold
% 	find(log.Measurement_Target==1 &  log.Measurement_NoiseDiode==0); % Hot
% 	find(log.Measurement_Target==5 &  log.Measurement_NoiseDiode==0); % Antenna
% 	find(log.Measurement_Target==2 &  log.Measurement_NoiseDiode==1); % Cold
% 	find(log.Measurement_Target==1 &  log.Measurement_NoiseDiode==1); % Hot
% 	find(log.Measurement_Target==5 &  log.Measurement_NoiseDiode==1); % Antenna
% 	};

% Checking the mean values to find the bug... From Axel
% Variations of mean amplitude and Tsys with time
TH=mean(log.T_Hot_Absorber);
drift.t  = log.t(hotIndices);
drift.T  = log.T_Hot_Absorber(hotIndices,:);
drift.a(1,:) = mean(rawSpectra(coldIndices,:),2,'omitnan');
drift.a(2,:) = mean(rawSpectra(hotIndices,:),2,'omitnan');
drift.a(3,:) = mean(rawSpectra(antennaIndices,:),2,'omitnan');
drift.Y    = drift.a(2,:) ./  drift.a(1,:);
drift.Tn   = (TH - drift.Y*TCold)./ (drift.Y-1);
drift.Ta   = (drift.a(3,:) - drift.a(1,:)) ./ (drift.a(2,:) - drift.a(1,:)) *(TH-TCold) + TCold;

figure
subplot(2,2,1); plot(drift.t, drift.Tn), ylabel('TN [K]') 
subplot(2,2,2); plot(drift.t, drift.Ta), ylabel('Ta [K]') 
subplot(2,2,3); plot(drift.t, drift.a),  ylabel('Counts [K]') 
subplot(2,2,4); plot(drift.t, drift.T),  ylabel('T Room  [K]') 
for i=1:4 subplot(2,2,i); set(gca, 'xlim', [0 24]); xlabel('t [h]'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Time interval for the calibration, splitting the indices for this day
% into the different calibration cycle



% Threshold for the separation, calibTime has to be in [min]
timeThresh=0:calibTime/60:24;

% Starting time for the complete cycle (will be the basis for separating
% the cycles according to calibrationTime). 
startingTimesHot=log.t(hotIndices);
startingTimesCold=log.t(coldIndices);
startingTimesAntennaAll=log.t(antennaIndices);
% Storing the indices specific to each calibration cycle in a new
% structure because by separating by time, we do not have the same
% number of individual cycle per calibration cycle
indices=struct();
for i = 1:length(timeThresh)-2
    condAntenna=startingTimesAntennaAll>timeThresh(i) & startingTimesAntennaAll<timeThresh(i+1);
    condHot=startingTimesHot>timeThresh(i) & startingTimesHot<timeThresh(i+1);
    condCold=startingTimesCold>timeThresh(i) & startingTimesCold<timeThresh(i+1);
   
    %indice=[validStartIndices(condHot); validStartIndices(condHot)+1; validStartIndices(condHot)+2; validStartIndices(condHot)+3; validStartIndices(condHot)+4; validStartIndices(condHot)+5];
    indices(i).validAntenna=antennaIndices(condAntenna);
    indices(i).validHot=hotIndices(condHot);
    indices(i).validCold=coldIndices(condCold);
end
condAntenna=startingTimesAntennaAll>timeThresh(length(timeThresh)-1);
condHot=startingTimesHot>timeThresh(length(timeThresh)-1);
condCold=startingTimesCold>timeThresh(length(timeThresh)-1);

%lastIndices=[validStartIndices(condHot); validStartIndices(condHot)+1; validStartIndices(condHot)+2; validStartIndices(condHot)+3; validStartIndices(condHot)+4; validStartIndices(condHot)+5];
%indices(length(timeThresh)-1).ind=lastIndices;
indices(length(timeThresh)-1).validAntenna=antennaIndices(condAntenna);
indices(length(timeThresh)-1).validHot=hotIndices(condHot);
indices(length(timeThresh)-1).validCold=coldIndices(condCold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hot and Cold spectra
% Outlier detection for hot and cold raw spectra as well as for their pointing
% angle
%
% Structure containing all the information about the calibrated spectra
calibratedSpectra=struct();

% Number of calibration cycles for this day:
nCalibrationCycles=length(indices);

for i=1:nCalibrationCycles
    ih=reshape(indices(i).validHot,[],1);
    ic=reshape(indices(i).validCold,[],1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Checking and removing any spurious angle for hot and cold
    hotAngleCheck=~(log.Elevation_Angle(ih)==retrievalTool.elevationAngleHot);
    coldAngleCheck=~(log.Elevation_Angle(ic)==retrievalTool.elevationAngleCold);
    
    if sum(hotAngleCheck>0)
        ih=ih(~hotAngleCheck);
        calibratedSpectra(i).hotAngleRemoved=sum(hotAngleCheck);
    else
        calibratedSpectra(i).hotAngleRemoved=0;
    end
    if sum(coldAngleCheck>0)
        ic=ic(~coldAngleCheck);
        calibratedSpectra(i).coldAngleRemoved=sum(coldAngleCheck);
    else
        calibratedSpectra(i).coldAngleRemoved=0;
    end
    
    initSizeHot=length(ih);
    initSizeCold=length(ic);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We compute the median of the calibration spectra as well as its std
    % deviation. We then remove the spectra which contains too many
    % channels that are beyond the (median +- n*stdDev)
    
    % "Standard" hot and cold raw counts for this cycle: (we could use other
    % standard like a daily mean or whatever
    medianRawCountsHot=nanmedian(rawSpectra(ih,:),1);
    medianRawCountsCold=nanmedian(rawSpectra(ic,:),1);
    
    medStdDevThreshHot=abs((rawSpectra(ih,:)-medianRawCountsHot))>retrievalTool.hotSpectraNumberOfStdDev*nanstd(rawSpectra(ih,:),1);
    medStdDevThreshCold=abs((rawSpectra(ic,:)-medianRawCountsCold))>retrievalTool.coldSpectraNumberOfStdDev*nanstd(rawSpectra(ic,:),1);
    
    ih=ih(sum(medStdDevThreshHot,2)<retrievalTool.threshNumRawSpectraHot);
    ic=ic(sum(medStdDevThreshCold,2)<retrievalTool.threshNumRawSpectraCold);
    
    calibratedSpectra(i).spuriousHotSpectra=initSizeHot-length(ih);
    calibratedSpectra(i).spuriousColdSpectra=initSizeCold-length(ic);  
    
    % Number of spectra left:
    calibratedSpectra(i).hotSpectraRemaining=length(ih);
    calibratedSpectra(i).coldSpectraRemaining=length(ic);
    
    % Saving clean hot and cold indices for this cycle
    calibratedSpectra(i).hotInd=ih;
    calibratedSpectra(i).coldInd=ic;
    
    % Final mean hot and cold raw counts for this cycle:
    calibratedSpectra(i).meanHotSpectra=nanmean(rawSpectra(ih,:),1);
    calibratedSpectra(i).meanColdSpectra=nanmean(rawSpectra(ic,:),1);
   
    calibratedSpectra(i).meanHotSpectra2=mean(rawSpectra(ih,:),2,'omitnan')';
    calibratedSpectra(i).meanColdSpectra2=mean(rawSpectra(ic,:),2,'omitnan')';
    
    calibratedSpectra(i).globalYFactor=calibratedSpectra(i).meanHotSpectra2/calibratedSpectra(i).meanColdSpectra2;
    
    % Final std dev hot and cold raw counts for this cycle:
    calibratedSpectra(i).stdHotSpectra=nanstd(rawSpectra(ih,:),1);
    calibratedSpectra(i).stdColdSpectra=nanstd(rawSpectra(ic,:),1);
    
    % Hot temperature corresponding to the hot spectra (to all spectra ?)
    calibratedSpectra(i).THot=nanmean(log.T_Hot_Absorber(ih));
    calibratedSpectra(i).stdTHot=nanstd(log.T_Hot_Absorber(ih));
    
    % === T_sys calculation : ===
    calibratedSpectra(i).Yfactor = (calibratedSpectra(i).meanHotSpectra)./( calibratedSpectra(i).meanColdSpectra);
    
    calibratedSpectra(i).TSys    = (calibratedSpectra(i).THot - calibratedSpectra(i).Yfactor*TCold)./(calibratedSpectra(i).Yfactor-1);
    
    % mean relative difference
    %threshHot=retrievalTool.threshRawSpectraHot*ones(size(ih));
    %threshCold=retrievalTool.threshRawSpectraCold*ones(size(ic));
    
    %mrdHot=(rawSpectra(ih,:)-meanRawCountsHot)./meanRawCountsHot > threshHot;
    %mrdCold=(rawSpectra(ic,:)-meanRawCountsCold)./meanRawCountsCold > threshCold;
    
    %plot(meanRawCountsHot)
    %hold on
    %plot(rawSpectra(2406,:))
    %plot((rawSpectra(2407,:)-meanRawCountsHot)./meanRawCountsHot)
    %ylim([lowerLim,upperLim])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
switch calType
    case 'standard'
        % Based on the discussion with Axel, changes in the calibration:
        % Averaging hot and cold FFTS counts on the time interval. 
        for i=1:nCalibrationCycles
            calibratedSpectra(i).calibrationVersion=calibVersion;
            calibratedSpectra(i).calibrationTime=calibTime;
            
            % All antenna measurements
            ia=reshape(indices(i).validAntenna,[],1);
            
            % Checking for NaN in the antenna spectra and keeping only complete
            % spectra for the calibration:
            ia=ia(sum(isnan(rawSpectra(ia,:)),2)<1);
            
            % Saving the indices for the Antenna
            calibratedSpectra(i).antennaInd=ia;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Flaging and removing (from the mean spectra only) bad angles for the antenna
            antennaAngleCheck=abs(log.Elevation_Angle(calibratedSpectra(i).antennaInd)-retrievalTool.elevationAngleAntenna)>retrievalTool.elevationAngleTolerance;
            if any(antennaAngleCheck)==1
                calibratedSpectra(i).antennaAngleRemoved=sum(antennaAngleCheck);
            else
                calibratedSpectra(i).antennaAngleRemoved=0;
            end
            calibratedSpectra(i).antennaIndCleanAngle=calibratedSpectra(i).antennaInd(~antennaAngleCheck);
            calibratedSpectra(i).numberOfCleanAntennaAngle=length(calibratedSpectra(i).antennaIndCleanAngle);
            
             calibratedSpectra(i).globalAntennaCounts=mean(rawSpectra(antennaIndices,:),2,'omitnan');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType='standard';
            
            % Mean Antenna counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);

            % Calibration
            calibratedSpectra(i).Tb = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra); 
        end
end
end

