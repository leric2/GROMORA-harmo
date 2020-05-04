function [drift,calibratedSpectra] = calibrate_generic(rawSpectra,standardLog,calibrationTool,calType)
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
% ARGUMENTS     | INPUTS: - rawSpectra: Matrix of raw data (1 line per
%               |           cycle, channels as columns.
%               |         - standardLog: standardized log file
%               |         - calibrationTool
%               |         - calType: Calibration type to do
%               |               1. standard: mean antenna vs mean hot/cold 
%               |                  spectra
%               |               2. debug: perform the standard one and:
%               |                  - Mean Up/Down: mean AT up, mean AT Down
%               |                  vs mean hot/cold
%               |                  - all cycles, mean hc: all individual 
%               |                  cycle vs mean hot/cold, uncleaned for 
%               |                  pointing angle problems.
%               |
%               | OUTPUTS: - drift: structure containing mean values of all
%               |            channels.
%               |          - calibratedSpectra: "structure array" of
%               |            calibrated data
%               |
% CALLS         | find_up_down_cycle()
%               |
%               |
%               |
%==========================================================================

% Calibration version
calibVersion='1.0.0';

% CalibrationTime in Minute
calibTime=calibrationTool.calibrationTime;

initialIndices={
    find(standardLog.Position==calibrationTool.indiceHot & standardLog.Tipping_Curve_active==0);      % Hot
    find(standardLog.Position==calibrationTool.indiceAntenna & standardLog.Tipping_Curve_active==0);  % Antenna
    find(standardLog.Position==calibrationTool.indiceCold & standardLog.Tipping_Curve_active==0);     % Cold
};

% Now we read complete half cycle Up (c-a-h) and Down (h-a-c) separately
[firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(standardLog,calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building the drift structure
% --> mainly for visualization
% 
% From Axel
% Variations of mean amplitude and Tsys with time
% need equal number of cycles !

for i=1:length(initialIndices); 	cn(i) = length(initialIndices{i}); end;
for i=1:length(initialIndices); 	initialIndices{i}=initialIndices{i}(1:min(cn)); end;

dailyMeanTHot=mean(standardLog.T_Hot_Absorber);

drift.t  = standardLog.t(initialIndices{1});
drift.T  = standardLog.T_Hot_Absorber(initialIndices{1});
drift.a(1,:) = mean(rawSpectra(initialIndices{1},:),2,'omitnan');
drift.a(2,:) = mean(rawSpectra(initialIndices{2},:),2,'omitnan');
drift.a(3,:) = mean(rawSpectra(initialIndices{3},:),2,'omitnan');
drift.Y    = drift.a(1,:) ./  drift.a(3,:);
drift.Tn   = (dailyMeanTHot - drift.Y*calibrationTool.TCold)./ (drift.Y-1);
drift.TSysLog=standardLog.FE_T_Sys(initialIndices{1});
drift.Ta   = (drift.a(2,:) - drift.a(3,:)) ./ (drift.a(1,:) - drift.a(3,:)) *(dailyMeanTHot-calibrationTool.TCold) + calibrationTool.TCold;

drift.dailyMedianHotSpectra=median(rawSpectra(initialIndices{1},:),1,'omitnan');
drift.dailyStdHotSpectra=nanstd(rawSpectra(initialIndices{1},:));
drift.dailyMedianColdSpectra=median(rawSpectra(initialIndices{3},:),1,'omitnan');
drift.dailyStdColdSpectra=nanstd(rawSpectra(initialIndices{3},:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time interval for the calibration, splitting the indices for this day
% into the different calibration cycle

% Threshold for the separation, calibTime has to be in [min]
dt=calibTime/60; % in [h]
timeThresh=0:dt:24;

% Starting time for the complete cycle (will be the basis for separating
% the cycles according to calibrationTime). 
startingTimesHot=standardLog.t(initialIndices{1});
startingTimesAntennaAll=standardLog.t(initialIndices{2});
startingTimesCold=standardLog.t(initialIndices{3});

startingTimesUp=standardLog.t(firstIndHalfUp);
startingTimesDown=standardLog.t(firstIndHalfDown);

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
    
    % The same for Up and Down calibration
    condUp=startingTimesUp>timeThresh(i) & startingTimesUp<timeThresh(i)+dt;
    condDown=startingTimesDown>timeThresh(i) & startingTimesDown<timeThresh(i)+dt;
    
    indices(i).validColdStartUp=[firstIndHalfUp(condUp); firstIndHalfUp(condUp)+1; firstIndHalfUp(condUp)+2];
    indices(i).validHotStartDown=[firstIndHalfDown(condDown); firstIndHalfDown(condDown)+1; firstIndHalfDown(condDown)+2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure containing all the information about the calibrated spectra
calibratedSpectra=struct();

% Hot and Cold spectra
% Outlier detection for hot and cold raw spectra as well as for their pointing
% angle

% Number of calibration cycles for this day
nCalibrationCycles=length(indices);

for i=1:nCalibrationCycles
    calibratedSpectra(i).theoreticalStartTime=timeThresh(i);
    
    % hot and cold indices for this calib cycle
    ih=reshape(indices(i).validHot,1,[]);
    ic=reshape(indices(i).validCold,1,[]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Checking and removing any spurious angle for hot and cold
    hotAngleCheck=~(standardLog.Elevation_Angle(ih)==calibrationTool.elevationAngleHot);
    coldAngleCheck=~(standardLog.Elevation_Angle(ic)==calibrationTool.elevationAngleCold);
    
    if sum(hotAngleCheck>0)
        ih=ih(~hotAngleCheck);
        %calibratedSpectra(i).hotAngleRemoved=sum(hotAngleCheck);
%     else
%         calibratedSpectra(i).hotAngleRemoved=0;
    end
    
    if sum(coldAngleCheck>0)
        ic=ic(~coldAngleCheck);
        %calibratedSpectra(i).coldAngleRemoved=sum(coldAngleCheck);
%     else
%         calibratedSpectra(i).coldAngleRemoved=0;
    end
    
    initSizeHot=length(ih);
    initSizeCold=length(ic);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO outlier detection for hot and cold spectra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % We compute the median of the calibration spectra as well as its std
    % deviation. We then remove the spectra which contains too many
    % channels that are beyond the (median +- n*stdDev)
          
    % Using the daily median and stddev to check for spurious hot and cold
    % spectra during this cycle:   
    medStdDevThreshHot=abs((rawSpectra(ih,:)-drift.dailyMedianHotSpectra))>3*drift.dailyStdHotSpectra;
    medStdDevThreshCold=abs((rawSpectra(ic,:)-drift.dailyMedianColdSpectra))>3*drift.dailyStdColdSpectra;
    
    %medStdDevThreshHot=abs((rawSpectra(ih,:)-drift.dailyMedianHotSpectra))>calibrationTool.hotSpectraNumberOfStdDev*drift.dailyStdHotSpectra;
    %medStdDevThreshCold=abs((rawSpectra(ic,:)-medianRawCountsCold))>calibrationTool.coldSpectraNumberOfStdDev*nanstd(rawSpectra(ic,:),1);
    
    ih=ih(sum(medStdDevThreshHot,2)<calibrationTool.threshNumRawSpectraHot);
    ic=ic(sum(medStdDevThreshCold,2)<calibrationTool.threshNumRawSpectraCold);
    
    calibratedSpectra(i).spuriousHotSpectra=initSizeHot-length(ih);
    calibratedSpectra(i).spuriousColdSpectra=initSizeCold-length(ic);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(ih)
        % Use drift structure for additionnal quality check
        k=find(drift.t >= standardLog.t(ih(1)) & drift.t <=standardLog.t(ih(end)));
    
        Tn_drift_i=drift.Tn(k);
        Ta_drift_i=drift.Ta(k);
    
        outliers = find( abs(Tn_drift_i-median(Tn_drift_i))>3*std(Tn_drift_i) | abs(Ta_drift_i-median(Ta_drift_i))>6*std(Ta_drift_i));
    
        if ~isempty(outliers)
%             out=[];
%             for n=1:length(outliers)
%                 out=[out ih(standardLog.t(ih) == drift.t(k(outliers(n))))];
%             end
            ih(outliers)=[];
            k(outliers)=[];
        end
        % Also used for stddev TSYS
        calibratedSpectra(i).TSysDrift=drift.Tn(k);
    else
        calibratedSpectra(i).TSysDrift=NaN;
    end
    calibratedSpectra(i).meanTSysDrift=nanmean(drift.Tn(k));
    calibratedSpectra(i).stdTSys=nanstd(drift.Tn(k));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO outlier detection for hot and cold spectra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number of spectra left:
    %calibratedSpectra(i).hotSpectraRemaining=length(ih);
    %calibratedSpectra(i).coldSpectraRemaining=length(ic);
    
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
    calibratedSpectra(i).THot=nanmean(standardLog.T_Hot_Absorber([ih,ic]));
    calibratedSpectra(i).stdTHot=nanstd(standardLog.T_Hot_Absorber([ih,ic]));
    
    % Computation of Final (clean) Tsys and its std deviation for this
    % cycle
    calibratedSpectra(i).Yspectral=calibratedSpectra(i).meanHotSpectra./calibratedSpectra(i).meanColdSpectra;
    
    % Noise Temperature Spectra (uncleaned for outliers in Y)
    calibratedSpectra(i).TN=(calibratedSpectra(i).THot - calibratedSpectra(i).Yspectral*calibrationTool.TCold)./ (calibratedSpectra(i).Yspectral -1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration of the antenna measurement with hot-cold load
switch calType
    case 'standard'
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
            
            % Saving all the indices for the Antenna
            calibratedSpectra(i).antennaInd=ia;
            calibratedSpectra(i).antennaIndUp=iaUp;
            calibratedSpectra(i).antennaIndDown=iaDown;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Flaging and removing (from the mean spectra only) bad angles for the antenna
            antennaAngleCheck=abs(standardLog.Elevation_Angle(calibratedSpectra(i).antennaInd)-calibrationTool.elevationAngleAntenna)>calibrationTool.elevationAngleTolerance;
            if sum(antennaAngleCheck)>0
                ia=ia(~antennaAngleCheck);
                %calibratedSpectra(i).antennaAngleRemoved=sum(antennaAngleCheck); 
            %else
                %calibratedSpectra(i).antennaAngleRemoved=0;
            end

            % Introduce here Outlier detection on antenna counts ?
            %TODO
            % ia=ia(no_outlier)
            
            calibratedSpectra(i).antennaIndCleanAngle=ia;
            %calibratedSpectra(i).numberOfCleanAntennaAngle=length(ia);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType=calType;
            
            % Antenna counts on this cycle:
            rsAntennaAll=rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:);
            
            % Mean Antenna counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);

            % Calibration
            % For every antenna measurement (to get the stddev of Tb on the
            % calibration cycle) against the mean of hot and cold spectra
            % for this calibration cycle.
            TbAll = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            calibratedSpectra(i).stdTb=std(TbAll);
            
            calibratedSpectra(i).Tb = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra); 
        end
    case 'debug'
        % In this mode, we compute 3 types of calibration for debugging
        % purpose:
        % 1. "standard" one: mean antenna vs mean hot/cold spectra
        % 2. "Mean Up/Down": mean AT up, mean AT Down vs mean hot/cold
        % 3. "all cycles, mean hc": all individual cycle vs mean hot/cold,
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
            antennaAngleCheck=abs(standardLog.Elevation_Angle(calibratedSpectra(i).antennaInd)-calibrationTool.elevationAngleAntenna)>calibrationTool.elevationAngleTolerance;
            if sum(antennaAngleCheck)>0
                ia=ia(~antennaAngleCheck);
                %calibratedSpectra(i).antennaAngleRemoved=sum(antennaAngleCheck); 
            %else
                %calibratedSpectra(i).antennaAngleRemoved=0;
            end

            % Introduce here Outlier detection on antenna counts ?
            %TODO
            % ia=ia(no_outlier)
            
            calibratedSpectra(i).antennaIndCleanAngle=ia;
            %calibratedSpectra(i).numberOfCleanAntennaAngle=length(ia);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Doing the calibration globally for this calibration cycle:
            calibratedSpectra(i).calibrationType=calType;
            
            % All clean Antenna counts on this cycle to compute the std Tb
            rsAntennaAllClean=rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:);
            
            TbAllCleanAngle = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaAllClean-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            calibratedSpectra(i).stdTb=std(TbAllCleanAngle);
            
            % To check the difference if we do the averaging after the
            % calibration of individual cycle:
            calibratedSpectra(i).meanFromTbAll=nanmean(TbAllCleanAngle);
            
            % Mean Antenna counts for this cycle
            rsAntenna=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanAngle,:),1);
            
            % Calibration 1. "standard"
            calibratedSpectra(i).Tb = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntenna-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra); 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing the calibration separately (Mean Up and Down) for this calibration cycle:
            angleCheckUp=abs(standardLog.Elevation_Angle(calibratedSpectra(i).antennaIndUp)-calibrationTool.elevationAngleAntenna)<calibrationTool.elevationAngleTolerance;
            angleCheckDown=abs(standardLog.Elevation_Angle(calibratedSpectra(i).antennaIndDown)-calibrationTool.elevationAngleAntenna)<calibrationTool.elevationAngleTolerance;
            
            % Cleaning the angle for "Mean Up/Down" calibration
            if any(angleCheckUp)==1
                calibratedSpectra(i).antennaUpAngleRemoved=sum(angleCheckUp); 
            else
                calibratedSpectra(i).antennaUpAngleRemoved=0;
            end
            calibratedSpectra(i).antennaIndCleanUpAngle=calibratedSpectra(i).antennaIndUp(~angleCheckUp)';
            calibratedSpectra(i).numberOfCleanAntennaUpAngle=length(calibratedSpectra(i).antennaIndCleanUpAngle);
            
            % Cleaning the angle for "Mean Up/Down" calibration
            if any(angleCheckDown)==1
                calibratedSpectra(i).antennaDownAngleRemoved=sum(angleCheckDown); 
            else
                calibratedSpectra(i).antennaDownAngleRemoved=0;
            end
            calibratedSpectra(i).antennaIndCleanDownAngle=calibratedSpectra(i).antennaIndDown(~angleCheckDown)';
            calibratedSpectra(i).numberOfCleanAntennaDownAngle=length(calibratedSpectra(i).antennaIndCleanDownAngle);
            
            % Mean Antenna Up/Down for this cycle
            rsAntennaUp=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanUpAngle,:),1);
            rsAntennaDown=nanmean(rawSpectra(calibratedSpectra(i).antennaIndCleanDownAngle,:),1);
            
            % Calibration 2. "Mean Up/Down"
            calibratedSpectra(i).TbUp = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaUp-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            calibratedSpectra(i).TbDown = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaDown-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calibration of all cycle separately. 
            % Saving the cycle by cycle calibration (separated for clarity)
            
            % For now, we do not keep the difference between up/dow cycle.
            % If one exist, it should be visible from "Mean Up/Down" type
            % of calibration. 
            rsAntennaAll=rawSpectra(calibratedSpectra(i).antennaInd,:);
            %rsAntennaUpAll=rawSpectra(calibratedSpectra(i).antennaIndUp,:);
            %rsAntennaDownAll=rawSpectra(calibratedSpectra(i).antennaIndDown,:);
            
            % Calibration 3. "all cycles, mean hc"
            calibratedSpectra(i).TbAll = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            % Angle flag all
            calibratedSpectra(i).antennaAngleFlagAll=antennaAngleCheck';
            
            %calibratedSpectra(i).TbUpAll = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaUpAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            % Angle flag up
            %calibratedSpectra(i).antennaAngleFlagUp=angleCheckUp';
            
            %calibratedSpectra(i).TbDownAll = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntennaDownAll-calibratedSpectra(i).meanColdSpectra)./(calibratedSpectra(i).meanHotSpectra-calibratedSpectra(i).meanColdSpectra);
            % Angle flag down
            %calibratedSpectra(i).antennaAngleFlagDown=angleCheckDown';
            
            % For the future: propagate the initial uncertainties into the
            % calibration formula.
        end
%     case 'time'
%         nCalibrationCycles=length(indices);
%         % Based on the starting times, we will then group the cycles together to
%         % form x minutes calibrationTime
%         %
%         %rsAntenna=ones(nCalibrationCycles,calibrationTool.numberOfChannels)*NaN;
%         %rsCold=ones(nCalibrationCycles,calibrationTool.numberOfChannels)*NaN;
%         %rsHot=ones(nCalibrationCycles,calibrationTool.numberOfChannels)*NaN;
%         %THot=ones(nCalibrationCycles,1)*NaN;
%         %stdTHot=ones(nCalibrationCycles,1)*NaN;
%         
%         % We need to loop through the calibration cycles because the number of averaged
%         % spectra might be different between each calibration cycle.
%         for i=1:nCalibrationCycles
%             calibratedSpectra(i).calibrationVersion=calibVersion;
%             calibratedSpectra(i).startInd=indices(i).ind(1,1);
%             calibratedSpectra(i).calibrationTime=calibTime;
%             
%             ih=reshape(indices(i).ind([1,6],:),[],1);
%             ia=reshape(indices(i).ind([2,5],:),[],1);
%             ic=reshape(indices(i).ind([3,4],:),[],1);
%             
%             % Computing some useful quantities for this cycles:
%             % Mean hot temperature for this cycle as well as its standard
%             % deviation:
%             calibratedSpectra(i).THot=nanmean(standardLog.T_Hot_Absorber(reshape(indices(i).ind,[],1))); 
%             calibratedSpectra(i).stdTHot=nanstd(standardLog.T_Hot_Absorber(reshape(indices(i).ind,[],1)));
%             
%             % Number of hot/cold/antenna averaged spectra for this cycle
%             % Considering all spectra that are not 100% NaN ... 
%             calibratedSpectra(i).nAvgSpectraHot=length(ih)-sum(all(isnan(rawSpectra(ih,:)),1));
%             calibratedSpectra(i).nAvgSpectraAntenna=length(ia)-sum(all(isnan(rawSpectra(ia,:)),1));
%             calibratedSpectra(i).nAvgSpectraCold=length(ic)-sum(all(isnan(rawSpectra(ic,:)),1));
%             
%             % Here we store the complete list of indices that were
%             % considered in this calibration cycle. It enables us to
%             % retrieve all standardLog data later and make the quality checks in a
%             % dedicated functions.
%             
%             % We take only the spectra that are not 100% NaN...
%             if sum(all(isnan(rawSpectra(ih,:)),1))==0
%                 calibratedSpectra(i).hotInd=ih;
%             else
%                 % TODO
%             end
%             if sum(all(isnan(rawSpectra(ia,:)),1))==0
%                 calibratedSpectra(i).antennaInd=ia;
%             else
%                 % TODO
%             end
%             if sum(all(isnan(rawSpectra(ic,:)),1))==0
%                 calibratedSpectra(i).coldInd=ic;
%             else
%                 % TODO
%             end
%             
%             % Mean and stdDev Temperature of the system for this cycle
%             %calibratedSpectra(i).Tsys=nanmean(standardLog.FE_T_Sys(reshape(indices(i).ind,[],1)));
%             %calibratedSpectra(i).stdTSys=nanstd(standardLog.FE_T_Sys(reshape(indices(i).ind,[],1)));
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Doing the calibration for individual cycle inside this
%             % calibration cycle (with global Thot):
%             ihAll=[indices(i).ind(1,:);indices(i).ind(6,:)];
%             iaAll=[indices(i).ind(2,:);indices(i).ind(5,:)];
%             icAll=[indices(i).ind(3,:);indices(i).ind(3,:)];
%             
%             % Mean Antenna counts for individual cycle (nCalibrationCycles x #channels)
%             rsAntennaAll=nanmean(cat(3,rawSpectra(iaAll(1,:),:),rawSpectra(iaAll(2,:),:)),3);
%             
%             % Mean hot counts for individual cycle (nCalibrationCycles x #channels)
%             rsHotAll=nanmean(cat(3,rawSpectra(ihAll(1,:),:),rawSpectra(ihAll(2,:),:)),3);
%             
%             % Mean cold counts for each cycle (nCalibrationCycles x #channels)
%             rsColdAll=nanmean(cat(3,rawSpectra(icAll(1,:),:),rawSpectra(icAll(2,:),:)),3);
%         
%             % Calibration on individual cycle (not saved)
%             TbAll = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold)'.*(rsAntennaAll-rsColdAll)./(rsHotAll-rsColdAll);
%             stdTbAll =nanstd(TbAll);
%             
%             % Checks
%             % calibratedSpectra(i).meanTbAll=nanmean(TbAll);
% 
%             % TOCHECK, not saved yet
%             %calibratedSpectra(i).TbAll=TbAll;
%             %calibratedSpectra(i).stdTbAll=stdTbAll;
%                    
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Doing the calibration globally for this calibration cycle:
%             calibratedSpectra(i).calibrationType='Counts avg and then calibrated';
%             
%             % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
%             rsHot=nanmean(rawSpectra(ih,:),1);
%             rsAntenna=nanmean(rawSpectra(ia,:),1);
%             rsCold=nanmean(rawSpectra(ic,:),1);
% 
%             % Calibration
%             calibratedSpectra(i).Tb = calibrationTool.TCold + (calibratedSpectra(i).THot-calibrationTool.TCold).*(rsAntenna-rsCold)./(rsHot-rsCold);
%             %calibratedSpectra(i).stdTb=nanstd(calibratedSpectra(i).Tb);
%             
%         end
%     case 'all'
%         % TODO
%         % If no calibration time is provided, we calibrate every cycle
%         % (2-1-0-0-1-2) and we don't need to loop
%         indicesAll=[validStartIndices; validStartIndices+1; validStartIndices+2; validStartIndices+3; validStartIndices+4; validStartIndices+5];
%         % Number of calibration cycle for this day:
%         nCalibrationCycles=size(indicesAll,2);
%         
%         % Here we don't need a loop nor a struct to store our indices as
%         % their number is the same in every cycle (6).
%         % So we store it in 2D matrices:
%         ih=[indicesAll(1,:);indicesAll(6,:)];
%         ia=[indicesAll(2,:);indicesAll(5,:)];
%         ic=[indicesAll(3,:);indicesAll(4,:)];
%         
%         % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
%         rsAntenna=nanmean(cat(3,rawSpectra(ia(1,:),:),rawSpectra(ia(2,:),:)),3);
%         
%         % Mean hot counts for each cycle (nCalibrationCycles x #channels)
%         rsHot=nanmean(cat(3,rawSpectra(ih(1,:),:),rawSpectra(ih(2,:),:)),3);
%         dailyMeanTHot=nanmean(standardLog.T_Hot_Absorber(indicesAll),1);
%         stdTHot=nanstd(standardLog.T_Hot_Absorber(indicesAll),1);
%         
%         % Mean cold counts for each cycle (nCalibrationCycles x #channels)
%         rsCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
%         
%         % Calibration
%         Tb = calibrationTool.TCold + (dailyMeanTHot-calibrationTool.TCold)'.*(rsAntenna-rsCold)./(rsHot-rsCold);
%         
%         % Mean System Temperature for each cycle
%         %Tsys=nanmean(standardLog.FE_T_Sys(indices),1);
%         % Std deviation of System Temperature for each cycle
%         %stdTSys=nanstd(standardLog.FE_T_Sys(indices),1);
%         
%         calibratedSpectra=struct();
%         % And we fill the final structure for the calibrated spectra
%         for i=1:nCalibrationCycles
%             calibratedSpectra(i).calibrationType='all cycle calibrated';
%             calibratedSpectra(i).calibrationVersion=calibVersion;
%             calibratedSpectra(i).hotInd=ih(:,i);
%             calibratedSpectra(i).antennaInd=ia(:,i);
%             calibratedSpectra(i).coldInd=ic(:,i);
%             calibratedSpectra(i).Tb=Tb(i,:);
%             calibratedSpectra(i).THot=dailyMeanTHot(i);
%             calibratedSpectra(i).stdTHot=stdTHot(i);
%         end
%     case 'all_then_avg'standardLog.Position==calibrationTool.indiceHot
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
%         THot=nanmean(standardLog.T_Hot_Absorber(indices),1);
%         stdTHot=nanstd(standardLog.T_Hot_Absorber(indices),1);
%         
%         % Mean cold counts for each cycle (nCalibrationCycles x #channels)
%         rsCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
%         
%         % Calibration
%         Tb = calibrationTool.TCold + (THot-calibrationTool.TCold)'.*(rsAntenna-rsCold)./(rsHot-rsCold);
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

%     % Nested function to extract all completed cycle from a given day based
%     % on the standardLog structure of this day. 
%     function firstIndCompleteCycle = find_completed_cycle(standardLog,calibrationTool)
%         indCold=calibrationTool.indiceCold;
%         indHot=calibrationTool.indiceHot;
%         indAntenna=calibrationTool.indiceAntenna;
%         hotInd=find(standardLog.Tipping_Curve_active==0 & (standardLog.Position==indHot));
%         extendedPos=[standardLog.Position -9999 -9999 -9999 -9999 -9999];
%         firstIndCompleteCycle=hotInd((extendedPos(hotInd+1)==indAntenna & extendedPos(hotInd+2)==indCold & extendedPos(hotInd+3)==indCold & extendedPos(hotInd+4)==indAntenna & extendedPos(hotInd+5)==indHot));
%     end
% 
%     % UNUSED
%     function indices = find_indice(standardLog,type,calibrationTool)
%         switch type
%             case 'cold'
%                 ind=calibrationTool.indiceCold;
%             case 'hot'
%                 ind=calibrationTool.indiceHot;
%             case 'antenna' 
%                 ind=calibrationTool.indiceAntenna;
%             otherwise
%             error('No valid type of indices provided')
%         end
%         indices=find(standardLog.Position==ind & standardLog.Tipping_Curve_active==0);
%     end   

end
