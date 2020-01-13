function calibratedSpectra = calibrate_generic(rawSpectra,log,retrievalTool,TCold,calType)
%CALIBRATE_GENERIC Summary of this function goes herecalibratedSpectra
%   Detailed explanation goes here

% Calibration version
calibVersion='1.0.0';

% CalibrationTime in Minute
calibTime=retrievalTool.calibrationTime;

% We define a complete cycle as 2-1-0-0-1-2:
hotIndCompleteCycle=find_completed_cycle(log,retrievalTool);

% Now we read cycle Up and Down separately:
[firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(log,retrievalTool);

% Group all indices by type for all the completed cycles:
ih=reshape([hotIndCompleteCycle; hotIndCompleteCycle+5],[],1);
ia=reshape([hotIndCompleteCycle+1;hotIndCompleteCycle+4],[],1);
ic=reshape([hotIndCompleteCycle+2;hotIndCompleteCycle+3],[],1);

% Checking all angle position of this day to remove any cycle with wrong
% angle
flagAngle=0;
% Hot (== retrievalTool.elevationAngleHot)
if ~log.Elevation_Angle(ih) == retrievalTool.elevationAngleHot
    % Find the wrong one and delete remove it from the cycle ?
    warning('There is a wrong angle for the hot load')
    flagAngle=1;
end

% Antenna
if ~((log.Elevation_Angle(ia)>retrievalTool.elevationAngleAntenna-retrievalTool.elevationAngleTolerance) & (log.Elevation_Angle(ia)<retrievalTool.elevationAngleAntenna+retrievalTool.elevationAngleTolerance))
    % Find the wrong one and delete remove it from the cycle ?
    warning('There is a wrong angle for the Antenna') 
    flagAngle=1;
end

% Cold (== retrievalTool.elevationAngleCold)
if ~log.Elevation_Angle(ic) == retrievalTool.elevationAngleCold
    % Find the wrong one and delete remove it from the cycle ?
    warning('There is a wrong angle for the Cold load')
    flagAngle=1;
end

% Remove cycles with wrong angles  
if flagAngle==0
    validStartIndices=hotIndCompleteCycle;                  % TODO
    validStartUp=firstIndHalfUp;
    validStartDown=firstIndHalfDown;
end

% Starting time for the complete cycle (will be the basis for separating
% the cycles according to calibrationTime). 
startingTimes=log.t(validStartIndices);
startingTimesUp=log.t(validStartUp);
startingTimesDown=log.t(validStartDown);

timeThresh=0:calibTime/60:24;

% Storing the indices specific to each calibration cycle in a new
% structure because by separating by time, we do not have the same
% number of individual cycle per calibration cycle
indices=struct();
for i = 1:length(timeThresh)-2
    cond=startingTimes>timeThresh(i) & startingTimes<timeThresh(i+1);
    condUp=startingTimesUp>timeThresh(i) & startingTimesUp<timeThresh(i+1);
    condDown=startingTimesDown>timeThresh(i) & startingTimesDown<timeThresh(i+1);
    
    indice=[validStartIndices(cond); validStartIndices(cond)+1; validStartIndices(cond)+2; validStartIndices(cond)+3; validStartIndices(cond)+4; validStartIndices(cond)+5];
    indices(i).ind=indice;
    
    indices(i).validColdStartUp=[validStartUp(condUp); validStartUp(condUp)+1; validStartUp(condUp)+2];
    indices(i).validHotStartDown=[validStartDown(condDown); validStartDown(condDown)+1; validStartDown(condDown)+2];
end
cond=startingTimes>timeThresh(length(timeThresh)-1);
condUp=startingTimesUp>timeThresh(length(timeThresh)-1);
condDown=startingTimesDown>timeThresh(length(timeThresh)-1);

lastIndices=[validStartIndices(cond); validStartIndices(cond)+1; validStartIndices(cond)+2; validStartIndices(cond)+3; validStartIndices(cond)+4; validStartIndices(cond)+5];
indices(length(timeThresh)-1).ind=lastIndices;
indices(length(timeThresh)-1).validColdStartUp=[validStartUp(condUp); validStartUp(condUp)+1; validStartUp(condUp)+2];;
indices(length(timeThresh)-1).validHotStartDown=[validStartDown(condDown); validStartDown(condDown)+1; validStartDown(condDown)+2];

calibratedSpectra=struct();

switch calType
    case 'hot_cold_avg_all_antenna'
        % Based on the discussion with Axel, changes in the calibration:
        % Averaging hot and cold FFTS counts on the time interval. 
        % Keeping each individual cycle for later analysis but computing as
        % well the avg calibrated spectra. 
        
        % Introducing some quality checks on hot and cold spectra for
        % removing bad quality spectra --> based on std ?
        
        nCalibrationCycles=length(indices);
        
        for i=1:nCalibrationCycles
        
            ih=reshape(indices(i).ind([1,6],:),[],1);
            ia=reshape(indices(i).ind([2,5],:),[],1);
            ic=reshape(indices(i).ind([3,4],:),[],1);
            
            % mean hot and cold raw counts for this cycle:
            rsMeanCold=nanmean(rawSpectra(ic,:),1);
            rsMeanHot=nanmean(rawSpectra(ih,:),1);
            
            % mean relative difference ????????
            MRD=(rawSpectra(2406,:)-rsMeanHot)./rsMeanHot;
            
            plot(MRD)
            
            plot(rsMeanHot)
            plot(rawSpectra(2406,:))
            ylim([lowerLim,upperLim])
        
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
            calibratedSpectra(i).THot=nanmean(log.T_Hot(reshape(indices(i).ind,[],1))); 
            calibratedSpectra(i).stdTHot=nanstd(log.T_Hot(reshape(indices(i).ind,[],1)));
            
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
        THot=nanmean(log.T_Hot(indicesAll),1);
        stdTHot=nanstd(log.T_Hot(indicesAll),1);
        
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
    case 'all_then_avg'
        timeThresh=0:calibTime/60:24;
        
        % Storing the indices specific to each calibration cycle in a new
        % structure because by separating by time, we do not have the same
        % number of individual cycle per calibration cycle
        avgIndices=struct();
        for i = 1:length(timeThresh)-2
            cond=startingTimes>timeThresh(i) & startingTimes<timeThresh(i+1);
            indice=[validStartIndices(cond); validStartIndices(cond)+1; validStartIndices(cond)+2; validStartIndices(cond)+3; validStartIndices(cond)+4; validStartIndices(cond)+5];
            avgIndices(i).ind=indice;
        end
        cond=startingTimes>timeThresh(length(timeThresh)-1);
        lastIndices=[validStartIndices(cond); validStartIndices(cond)+1; validStartIndices(cond)+2; validStartIndices(cond)+3; validStartIndices(cond)+4; validStartIndices(cond)+5];
        avgIndices(length(timeThresh)-1).ind=lastIndices;
        
        nCalibrationCycles=length(avgIndices);
        
        indices=[validStartIndices; validStartIndices+1; validStartIndices+2; validStartIndices+3; validStartIndices+4; validStartIndices+5];

        % Here we don't need a loop nor a struct to store our indices as
        % their number is the same in every cycle (6).
        % So we store it in 2D matrices:
        ih=[indices(1,:);indices(6,:)];
        ia=[indices(2,:);indices(5,:)];
        ic=[indices(3,:);indices(4,:)];
        
        % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
        rsAntenna=nanmean(cat(3,rawSpectra(ia(1,:),:),rawSpectra(ia(2,:),:)),3);
        
        % Mean hot counts for each cycle (nCalibrationCycles x #channels)
        rsHot=nanmean(cat(3,rawSpectra(ih(1,:),:),rawSpectra(ih(2,:),:)),3);
        THot=nanmean(log.T_Hot(indices),1);
        stdTHot=nanstd(log.T_Hot(indices),1);
        
        % Mean cold counts for each cycle (nCalibrationCycles x #channels)
        rsCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
        
        % Calibration
        Tb = TCold + (THot-TCold)'.*(rsAntenna-rsCold)./(rsHot-rsCold);
        
        % Now doing the averaging of the brightness temperature
        % "calibTime"
        startIndiceCycle=1;
        for i = 1:nCalibrationCycles-1
            indiceForThisCycle=startIndiceCycle:startIndiceCycle+size(avgIndices(i).ind,2);
            %indiceForThisCycle=reshape(avgIndices(i).ind,[],1);
            calibratedSpectra(i).calibrationType='all cycle calibrated and then avg';
            calibratedSpectra(i).calibrationVersion=calibVersion;
            calibratedSpectra(i).hotInd=reshape(avgIndices(i).ind([1,6],:),[],1);
            calibratedSpectra(i).antennaInd=reshape(avgIndices(i).ind([2,5],:),[],1);
            calibratedSpectra(i).coldInd=reshape(avgIndices(i).ind([3,4],:),[],1);
            calibratedSpectra(i).Tb=nanmean(Tb(indiceForThisCycle,:),1);
            calibratedSpectra(i).THot=nanmean(THot(indiceForThisCycle),1);
            startIndiceCycle=startIndiceCycle+length(avgIndices(i).ind);
            calibratedSpectra(i).stdTHot=0;
        end
        indiceForThisCycle=startIndiceCycle:size(Tb,1);
        %indiceForThisCycle=reshape(avgIndices(i).ind,[],1);
        calibratedSpectra(48).calibrationType='all cycle calibrated and then avg';
        calibratedSpectra(48).calibrationVersion=calibVersion;
        calibratedSpectra(48).hotInd=reshape(avgIndices(i).ind([1,6],:),[],1);
        calibratedSpectra(48).antennaInd=reshape(avgIndices(i).ind([2,5],:),[],1);
        calibratedSpectra(48).coldInd=reshape(avgIndices(i).ind([3,4],:),[],1);
        calibratedSpectra(48).Tb=nanmean(Tb(indiceForThisCycle,:),1);
        calibratedSpectra(48).THot=nanmean(THot(indiceForThisCycle),1);
        calibratedSpectra(48).stdTHot=0;
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

