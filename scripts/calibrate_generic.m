function [drift,calibratedSpectra] = calibrate_generic(rawSpectra,log,retrievalTool,TCold,calType)
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
%hotIndices=find(log.Position==retrievalTool.indiceHot & log.Tipping_Curve_active==0);
%antennaIndices=find(log.Position==retrievalTool.indiceAntenna & log.Tipping_Curve_active==0);
%coldIndices=find(log.Position==retrievalTool.indiceCold & log.Tipping_Curve_active==0);

initialIndices={
    find(log.Position==retrievalTool.indiceHot & log.Tipping_Curve_active==0);      % Hot
    find(log.Position==retrievalTool.indiceAntenna & log.Tipping_Curve_active==0);  % Antenna
    find(log.Position==retrievalTool.indiceCold & log.Tipping_Curve_active==0);     % Cold
};

% Now we read complete half cycle Up (c-a-h) and Down (h-a-c) separately
[firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(log,retrievalTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the mean values to find the bug... From Axel
% Variations of mean amplitude and Tsys with time
% ensure equal number of cycles

for i=1:length(initialIndices); 	cn(i) = length(initialIndices{i}); end;
for i=1:length(initialIndices); 	initialIndices{i}=initialIndices{i}(1:min(cn)); end;

TH=mean(log.T_Hot_Absorber);
drift.t  = log.t(initialIndices{1});
drift.T  = log.T_Hot_Absorber(initialIndices{1});
drift.a(1,:) = mean(rawSpectra(initialIndices{1},:),2,'omitnan');
drift.a(2,:) = mean(rawSpectra(initialIndices{2},:),2,'omitnan');
drift.a(3,:) = mean(rawSpectra(initialIndices{3},:),2,'omitnan');
drift.Y    = drift.a(1,:) ./  drift.a(3,:);
drift.Tn   = (TH - drift.Y*TCold)./ (drift.Y-1);
drift.TSysLog=log.FE_T_Sys(initialIndices{1});
drift.Ta   = (drift.a(2,:) - drift.a(3,:)) ./ (drift.a(1,:) - drift.a(3,:)) *(TH-TCold) + TCold;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Time interval for the calibration, splitting the indices for this day
% into the different calibration cycle

% Threshold for the separation, calibTime has to be in [min]
dt=calibTime/60; % in [h]
timeThresh=0:dt:24;

% Starting time for the complete cycle (will be the basis for separating
% the cycles according to calibrationTime). 
startingTimesHot=log.t(initialIndices{1});
startingTimesAntennaAll=log.t(initialIndices{2});
startingTimesCold=log.t(initialIndices{3});

startingTimesUp=log.t(firstIndHalfUp);
startingTimesDown=log.t(firstIndHalfDown);

% Storing the indices specific to each calibration cycle in a new
% structure because by separating by time, we do not have the same
% number of individual cycle per calibration cycle
indices=struct();
for i = 1:length(timeThresh)-1
    %condAntenna=startingTimesAntennaAll>timeThresh(i) & startingTimesAntennaAll<timeThresh(i+1);
    condAntenna=startingTimesAntennaAll>timeThresh(i) & startingTimesAntennaAll<timeThresh(i)+dt;
    %condHot=startingTimesHot>timeThresh(i) & startingTimesHot<timeThresh(i+1);
    condHot=startingTimesHot>timeThresh(i) & startingTimesHot<timeThresh(i)+dt;
    %condCold=startingTimesCold>timeThresh(i) & startingTimesCold<timeThresh(i+1);
    condCold=startingTimesCold>timeThresh(i) & startingTimesCold<timeThresh(i)+dt;
    
    %condUp=startingTimesUp>timeThresh(i) & startingTimesUp<timeThresh(i+1);
    condUp=startingTimesUp>timeThresh(i) & startingTimesUp<timeThresh(i)+dt;
    %condDown=startingTimesDown>timeThresh(i) & startingTimesDown<timeThresh(i+1);
    condDown=startingTimesDown>timeThresh(i) & startingTimesDown<timeThresh(i)+dt;
    
    %indice=[validStartIndices(condHot); validStartIndices(condHot)+1; validStartIndices(condHot)+2; validStartIndices(condHot)+3; validStartIndices(condHot)+4; validStartIndices(condHot)+5];
    indices(i).validAntenna=initialIndices{2}(condAntenna);
    indices(i).validHot=initialIndices{1}(condHot);
    indices(i).validCold=initialIndices{3}(condCold);
    
    %indices(i).ind=indice;
    indices(i).validColdStartUp=[firstIndHalfUp(condUp); firstIndHalfUp(condUp)+1; firstIndHalfUp(condUp)+2];
    indices(i).validHotStartDown=[firstIndHalfDown(condDown); firstIndHalfDown(condDown)+1; firstIndHalfDown(condDown)+2];
end

% condAntenna=startingTimesAntennaAll>timeThresh(length(timeThresh)-1);
% condHot=startingTimesHot>timeThresh(length(timeThresh)-1);
% condCold=startingTimesCold>timeThresh(length(timeThresh)-1);
% condUp=startingTimesUp>timeThresh(length(timeThresh)-1);
% condDown=startingTimesDown>timeThresh(length(timeThresh)-1);

%lastIndices=[validStartIndices(condHot); validStartIndices(condHot)+1; validStartIndices(condHot)+2; validStartIndices(condHot)+3; validStartIndices(condHot)+4; validStartIndices(condHot)+5];
%indices(length(timeThresh)-1).ind=lastIndices;
% indices(length(timeThresh)-1).validAntenna=antennaIndices(condAntenna);
% indices(length(timeThresh)-1).validHot=hotIndices(condHot);
% indices(length(timeThresh)-1).validCold=coldIndices(condCold);
% indices(length(timeThresh)-1).validColdStartUp=[firstIndHalfUp(condUp); firstIndHalfUp(condUp)+1; firstIndHalfUp(condUp)+2];
% indices(length(timeThresh)-1).validHotStartDown=[firstIndHalfDown(condDown); firstIndHalfDown(condDown)+1; firstIndHalfDown(condDown)+2];

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
    calibratedSpectra(i).theoreticalStartTime=timeThresh(i);
    
    ih=reshape(indices(i).validHot,1,[]);
    ic=reshape(indices(i).validCold,1,[]);

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
    
    ihc=sort([ih,ic]);
    
    % Use drift structure for quality check ?
    k=find(drift.t>=log.t(ihc(1)) & drift.t<=log.t(ihc(1))+dt);
    
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
    
    % Final std dev hot and cold raw counts for this cycle:
    calibratedSpectra(i).stdHotSpectra=nanstd(rawSpectra(ih,:),1);
    calibratedSpectra(i).stdColdSpectra=nanstd(rawSpectra(ic,:),1);
    
    % Hot temperature corresponding to the hot spectra (to all spectra ?)
    calibratedSpectra(i).THot=nanmean(log.T_Hot_Absorber([ih,ic]));
    calibratedSpectra(i).stdTHot=nanstd(log.T_Hot_Absorber([ih,ic]));
    
    % Computation of Final (clean) Tsys and its std deviation for this cycle
%     if length(ih) == length(ic)
%         YAll=rawSpectra(ih,:)./rawSpectra(ic,:);
%         
%     else
%         
%     end
    %calibratedSpectra(i).Ymean=nanmean(calibratedSpectra(i).meanHotSpectra)/mean(calibratedSpectra(i).meanColdSpectra);
    
    %calibratedSpectra(i).TSys=(calibratedSpectra(i).THot-TCold*calibratedSpectra(i).Ymean)/(calibratedSpectra(i).Ymean-1);
    %calibratedSpectra(i).stdTSys=nanstd(calibratedSpectra(i).TSys);
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
            ia=reshape(indices(i).validAntenna,1,[]);
            
            % Antenna measurements inside a half cycle
            iaUp=reshape(indices(i).validColdStartUp(2,:),1,[]);
            iaDown=reshape(indices(i).validHotStartDown(2,:),1,[]);
            
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
            % Flaging and removing (from the mean spectra only) bad angles for the antenna
            antennaAngleCheck=abs(log.Elevation_Angle(calibratedSpectra(i).antennaInd)-retrievalTool.elevationAngleAntenna)>retrievalTool.elevationAngleTolerance;
            if any(antennaAngleCheck)==1
                calibratedSpectra(i).antennaAngleRemoved=sum(antennaAngleCheck);
            else
                calibratedSpectra(i).antennaAngleRemoved=0;
            end
            calibratedSpectra(i).antennaIndCleanAngle=calibratedSpectra(i).antennaInd(~antennaAngleCheck)';
            calibratedSpectra(i).numberOfCleanAntennaAngle=length(calibratedSpectra(i).antennaIndCleanAngle);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType='standard';
            
            % Mean Antenna counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);
            
            % Y factor
            calibratedSpectra(i).Y=calibratedSpectra(i).meanHotSpectra ./ calibratedSpectra(i).meanColdSpectra;
            
            % Noise Temperature Spectra
            calibratedSpectra(i).TN=(calibratedSpectra(i).THot - calibratedSpectra(i).Y*TCold)./ (calibratedSpectra(i).Y -1);
            
            % Calibration
            calibratedSpectra(i).Tb = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra); 
        end
    case 'debug'
        % Based on the discussion with Axel, changes in the calibration:
        % Averaging hot and cold FFTS counts on the time interval.
        % Keeping each individual cycle for later analysis but computing as
        % well the avg calibrated spectra.
        
        for i=1:nCalibrationCycles
            calibratedSpectra(i).calibrationVersion=calibVersion;
            calibratedSpectra(i).calibrationTime=calibTime;
            
            % All antenna measurements
            ia=reshape(indices(i).validAntenna,[],1);
            
            % Antenna measurements inside a half cycle
            iaUp=reshape(indices(i).validColdStartUp(2,:),[],1);
            iaDown=reshape(indices(i).validHotStartDown(2,:),[],1);
            
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
            % Flaging and removing (from the mean spectra only) bad angles for the antenna
            antennaAngleCheck=abs(log.Elevation_Angle(calibratedSpectra(i).antennaInd)-retrievalTool.elevationAngleAntenna)>retrievalTool.elevationAngleTolerance;
            if any(antennaAngleCheck)==1
                calibratedSpectra(i).antennaAngleRemoved=sum(antennaAngleCheck);
            else
                calibratedSpectra(i).antennaAngleRemoved=0;
            end
            calibratedSpectra(i).antennaIndCleanAngle=calibratedSpectra(i).antennaInd(~antennaAngleCheck);
            calibratedSpectra(i).numberOfCleanAntennaAngle=length(calibratedSpectra(i).antennaIndCleanAngle);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType='debug';
            
            % Mean Antenna counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);
            
            % Calibration
            calibratedSpectra(i).Tb = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration separately (Up and Down) for this calibration cycle:
            angleCheckUp=abs(log.Elevation_Angle(calibratedSpectra(i).antennaIndUp)-retrievalTool.elevationAngleAntenna)<retrievalTool.elevationAngleTolerance;
            angleCheckDown=abs(log.Elevation_Angle(calibratedSpectra(i).antennaIndDown)-retrievalTool.elevationAngleAntenna)<retrievalTool.elevationAngleTolerance;
            
            % Antenna Up for this cycle
            rsAntennaUp=nanmean(rawSpectra(calibratedSpectra(i).antennaIndUp,:),1);
            rsAntennaDown=nanmean(rawSpectra(calibratedSpectra(i).antennaIndDown,:),1);
            
            calibratedSpectra(i).TbUp = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntennaUp-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            calibratedSpectra(i).TbDown = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntennaDown-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Saving the cycle by cycle calibration (separated for clarity)
            rsAntennaAll=rawSpectra(calibratedSpectra(i).antennaInd,:);
            rsAntennaUpAll=rawSpectra(calibratedSpectra(i).antennaIndUp,:);
            rsAntennaDownAll=rawSpectra(calibratedSpectra(i).antennaIndDown,:);
            
            calibratedSpectra(i).TbAll = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntennaAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            % Angle flag
            % calibratedSpectra(i).antennaAngleFlagAll=angleCheck';
            
            calibratedSpectra(i).TbUpAll = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntennaUpAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            % Angle flag
            calibratedSpectra(i).antennaAngleFlagUp=angleCheckUp';
            
            calibratedSpectra(i).TbDownAll = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntennaDownAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            % Angle flag
            calibratedSpectra(i).antennaAngleFlagDown=angleCheckDown';
            
            % And the standard deviation of final Tb
            calibratedSpectra(i).stdTb=nanstd(calibratedSpectra(i).TbAll);
            
            calibratedSpectra(i).meanFromTbAll=mean(calibratedSpectra(i).TbAll);
            calibratedSpectra(i).meanFromTbUpAll=mean(calibratedSpectra(i).TbUpAll);
            calibratedSpectra(i).meanFromTbDownAll=mean(calibratedSpectra(i).TbDownAll);
            % For the future: propagate the initial uncertainties into the
            % calibration formula.
        end
    case 'time'
        nCalibrationCycles=length(indices);
        % Based on the starting times, we will then group the cycles together to
        % form x minutes calibrationTime
        %
        %rsAntenna=ones(nCalibrationCycles,retrievalTool.numberOfChannels)*NaN;
        %rsCold=ones(nCalibrationCycles,retrievalTool.numberOfChannels)*NaN;
        %rsHot=ones(nCalibrationCycles,retrievalTool.numberOfChannels)*NaN;
        %THot=ones(nCalibrationCycles,1)*NaN;
        %stdTHot=ones(nCalibrationCycles,1)*NaN;
        
        % We need to loop through the calibration cycles because the number of averaged
        % spectra might be different between each calibration cycle.
        for i=1:nCalibrationCycles
            calibratedSpectra(i).calibrationVersion=calibVersion;
            calibratedSpectra(i).startInd=indices(i).ind(1,1);
            calibratedSpectra(i).calibrationTime=calibTime;
            
            ih=reshape(indices(i).ind([1,6],:),[],1);
            ia=reshape(indices(i).ind([2,5],:),[],1);
            ic=reshape(indices(i).ind([3,4],:),[],1);
            
            % Computing some useful quantities for this cycles:
            % Mean hot temperature for this cycle as well as its standard
            % deviation:
            calibratedSpectra(i).THot=nanmean(log.T_Hot_Absorber(reshape(indices(i).ind,[],1))); 
            calibratedSpectra(i).stdTHot=nanstd(log.T_Hot_Absorber(reshape(indices(i).ind,[],1)));
            
            % Number of hot/cold/antenna averaged spectra for this cycle
            % Considering all spectra that are not 100% NaN ... 
            calibratedSpectra(i).nAvgSpectraHot=length(ih)-sum(all(isnan(rawSpectra(ih,:)),1));
            calibratedSpectra(i).nAvgSpectraAntenna=length(ia)-sum(all(isnan(rawSpectra(ia,:)),1));
            calibratedSpectra(i).nAvgSpectraCold=length(ic)-sum(all(isnan(rawSpectra(ic,:)),1));
            
            % Here we store the complete list of indices that were
            % considered in this calibration cycle. It enables us to
            % retrieve all log data later and make the quality checks in a
            % dedicated functions.
            
            % We take only the spectra that are not 100% NaN...
            if sum(all(isnan(rawSpectra(ih,:)),1))==0
                calibratedSpectra(i).hotInd=ih;
            else
                % TODO
            end
            if sum(all(isnan(rawSpectra(ia,:)),1))==0
                calibratedSpectra(i).antennaInd=ia;
            else
                % TODO
            end
            if sum(all(isnan(rawSpectra(ic,:)),1))==0
                calibratedSpectra(i).coldInd=ic;
            else
                % TODO
            end
            
            % Mean and stdDev Temperature of the system for this cycle
            %calibratedSpectra(i).Tsys=nanmean(log.FE_T_Sys(reshape(indices(i).ind,[],1)));
            %calibratedSpectra(i).stdTSys=nanstd(log.FE_T_Sys(reshape(indices(i).ind,[],1)));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration for individual cycle inside this
            % calibration cycle (with global Thot):
            ihAll=[indices(i).ind(1,:);indices(i).ind(6,:)];
            iaAll=[indices(i).ind(2,:);indices(i).ind(5,:)];
            icAll=[indices(i).ind(3,:);indices(i).ind(3,:)];
            
            % Mean Antenna counts for individual cycle (nCalibrationCycles x #channels)
            rsAntennaAll=nanmean(cat(3,rawSpectra(iaAll(1,:),:),rawSpectra(iaAll(2,:),:)),3);
            
            % Mean hot counts for individual cycle (nCalibrationCycles x #channels)
            rsHotAll=nanmean(cat(3,rawSpectra(ihAll(1,:),:),rawSpectra(ihAll(2,:),:)),3);
            
            % Mean cold counts for each cycle (nCalibrationCycles x #channels)
            rsColdAll=nanmean(cat(3,rawSpectra(icAll(1,:),:),rawSpectra(icAll(2,:),:)),3);
        
            % Calibration on individual cycle (not saved)
            TbAll = TCold + (calibratedSpectra(i).THot-TCold)'.*(rsAntennaAll-rsColdAll)./(rsHotAll-rsColdAll);
            stdTbAll =nanstd(TbAll);
            
            % Checks
            % calibratedSpectra(i).meanTbAll=nanmean(TbAll);

            % TOCHECK, not saved yet
            %calibratedSpectra(i).TbAll=TbAll;
            %calibratedSpectra(i).stdTbAll=stdTbAll;
                   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType='Counts avg and then calibrated';
            
            % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
            rsHot=nanmean(rawSpectra(ih,:),1);
            rsAntenna=nanmean(rawSpectra(ia,:),1);
            rsCold=nanmean(rawSpectra(ic,:),1);

            % Calibration
            calibratedSpectra(i).Tb = TCold + (calibratedSpectra(i).THot-TCold).*(rsAntenna-rsCold)./(rsHot-rsCold);
            %calibratedSpectra(i).stdTb=nanstd(calibratedSpectra(i).Tb);
            
        end
    case 'all'
        % TODO
        % If no calibration time is provided, we calibrate every cycle
        % (2-1-0-0-1-2) and we don't need to loop
        indicesAll=[validStartIndices; validStartIndices+1; validStartIndices+2; validStartIndices+3; validStartIndices+4; validStartIndices+5];
        % Number of calibration cycle for this day:
        nCalibrationCycles=size(indicesAll,2);
        
        % Here we don't need a loop nor a struct to store our indices as
        % their number is the same in every cycle (6).
        % So we store it in 2D matrices:
        ih=[indicesAll(1,:);indicesAll(6,:)];
        ia=[indicesAll(2,:);indicesAll(5,:)];
        ic=[indicesAll(3,:);indicesAll(4,:)];
        
        % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
        rsAntenna=nanmean(cat(3,rawSpectra(ia(1,:),:),rawSpectra(ia(2,:),:)),3);
        
        % Mean hot counts for each cycle (nCalibrationCycles x #channels)
        rsHot=nanmean(cat(3,rawSpectra(ih(1,:),:),rawSpectra(ih(2,:),:)),3);
        THot=nanmean(log.T_Hot_Absorber(indicesAll),1);
        stdTHot=nanstd(log.T_Hot_Absorber(indicesAll),1);
        
        % Mean cold counts for each cycle (nCalibrationCycles x #channels)
        rsCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
        
        % Calibration
        Tb = TCold + (THot-TCold)'.*(rsAntenna-rsCold)./(rsHot-rsCold);
        
        % Mean System Temperature for each cycle
        %Tsys=nanmean(log.FE_T_Sys(indices),1);
        % Std deviation of System Temperature for each cycle
        %stdTSys=nanstd(log.FE_T_Sys(indices),1);
        
        calibratedSpectra=struct();
        % And we fill the final structure for the calibrated spectra
        for i=1:nCalibrationCycles
            calibratedSpectra(i).calibrationType='all cycle calibrated';
            calibratedSpectra(i).calibrationVersion=calibVersion;
            calibratedSpectra(i).hotInd=ih(:,i);
            calibratedSpectra(i).antennaInd=ia(:,i);
            calibratedSpectra(i).coldInd=ic(:,i);
            calibratedSpectra(i).Tb=Tb(i,:);
            calibratedSpectra(i).THot=THot(i);
            calibratedSpectra(i).stdTHot=stdTHot(i);
        end
%     case 'all_then_avg'log.Position==retrievalTool.indiceHot
%         timeThresh=0:calibTime/60:24;
%         
%         % Storing the indices specific to each calibration cycle in a new
%         % structure because by separating by time, we do not have the same
%         % number of individual cycle per calibration cycle
%         avgIndices=struct();
%         for i = 1:length(timeThresh)-2
%             condHot=startingTimesHot>timeThresh(i) & startingTimesHot<timeThresh(i+1);
%             indice=[validStartIndices(condHot); validStartIndices(condHot)+1; validStartIndices(condHot)+2; validStartIndices(condHot)+3; validStartIndices(condHot)+4; validStartIndices(condHot)+5];
%             avgIndices(i).ind=indice;
%         end
%         condHot=startingTimesHot>timeThresh(length(timeThresh)-1);
%         lastIndices=[validStartIndices(condHot); validStartIndices(condHot)+1; validStartIndices(condHot)+2; validStartIndices(condHot)+3; validStartIndices(condHot)+4; validStartIndices(condHot)+5];
%         avgIndices(length(timeThresh)-1).ind=lastIndices;
%         
%         nCalibrationCycles=length(avgIndices);
%         
%         indices=[validStartIndices; validStartIndices+1; validStartIndices+2; validStartIndices+3; validStartIndices+4; validStartIndices+5];
% 
%         % Here we don't need a loop nor a struct to store our indices as
%         % their number is the same in every cycle (6).
%         % So we store it in 2D matrices:
%         ih=[indices(1,:);indices(6,:)];
%         ia=[indices(2,:);indices(5,:)];
%         ic=[indices(3,:);indices(4,:)];
%         
%         % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
%         rsAntenna=nanmean(cat(3,rawSpectra(ia(1,:),:),rawSpectra(ia(2,:),:)),3);
%         
%         % Mean hot counts for each cycle (nCalibrationCycles x #channels)
%         rsHot=nanmean(cat(3,rawSpectra(ih(1,:),:),rawSpectra(ih(2,:),:)),3);
%         THot=nanmean(log.T_Hot_Absorber(indices),1);
%         stdTHot=nanstd(log.T_Hot_Absorber(indices),1);
%         
%         % Mean cold counts for each cycle (nCalibrationCycles x #channels)
%         rsCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
%         
%         % Calibration
%         Tb = TCold + (THot-TCold)'.*(rsAntenna-rsCold)./(rsHot-rsCold);
%         
%         % Now doing the averaging of the brightness temperature
%         % "calibTime"
%         startIndiceCycle=1;
%         for i = 1:nCalibrationCycles-1
%             indiceForThisCycle=startIndiceCycle:startIndiceCycle+size(avgIndices(i).ind,2);
%             %indiceForThisCycle=reshape(avgIndices(i).ind,[],1);
%             calibratedSpectra(i).calibrationType='all cycle calibrated and then avg';
%             calibratedSpectra(i).calibrationVersion=calibVersion;
%             calibratedSpectra(i).hotInd=reshape(avgIndices(i).ind([1,6],:),[],1);
%             calibratedSpectra(i).antennaInd=reshape(avgIndices(i).ind([2,5],:),[],1);
%             calibratedSpectra(i).coldInd=reshape(avgIndices(i).ind([3,4],:),[],1);
%             calibratedSpectra(i).Tb=nanmean(Tb(indiceForThisCycle,:),1);
%             calibratedSpectra(i).THot=nanmean(THot(indiceForThisCycle),1);
%             startIndiceCycle=startIndiceCycle+length(avgIndices(i).ind);
%             calibratedSpectra(i).stdTHot=0;
%         end
%         indiceForThisCycle=startIndiceCycle:size(Tb,1);
%         %indiceForThisCycle=reshape(avgIndices(i).ind,[],1);
%         calibratedSpectra(48).calibrationType='all cycle calibrated and then avg';
%         calibratedSpectra(48).calibrationVersion=calibVersion;
%         calibratedSpectra(48).hotInd=reshape(avgIndices(i).ind([1,6],:),[],1);
%         calibratedSpectra(48).antennaInd=reshape(avgIndices(i).ind([2,5],:),[],1);
%         calibratedSpectra(48).coldInd=reshape(avgIndices(i).ind([3,4],:),[],1);
%         calibratedSpectra(48).Tb=nanmean(Tb(indiceForThisCycle,:),1);
%         calibratedSpectra(48).THot=nanmean(THot(indiceForThisCycle),1);
%         calibratedSpectra(48).stdTHot=0;
end

    % Nested function to extract all completed cycle from a given day based
    % on the log structure of this day. 
    function firstIndCompleteCycle = find_completed_cycle(log,retrievalTool)
        indCold=retrievalTool.indiceCold;
        indHot=retrievalTool.indiceHot;
        indAntenna=retrievalTool.indiceAntenna;
        hotInd=find(log.Tipping_Curve_active==0 & (log.Position==indHot));
        extendedPos=[log.Position -9999 -9999 -9999 -9999 -9999];
        firstIndCompleteCycle=hotInd((extendedPos(hotInd+1)==indAntenna & extendedPos(hotInd+2)==indCold & extendedPos(hotInd+3)==indCold & extendedPos(hotInd+4)==indAntenna & extendedPos(hotInd+5)==indHot));
    end

    % UNUSED
    function indices = find_indice(log,type,retrievalTool)
        switch type
            case 'cold'
                ind=retrievalTool.indiceCold;
            case 'hot'
                ind=retrievalTool.indiceHot;
            case 'antenna' 
                ind=retrievalTool.indiceAntenna;
            otherwise
            error('No valid type of indices provided')
        end
        indices=find(log.Position==ind & log.Tipping_Curve_active==0);
    end   

end

