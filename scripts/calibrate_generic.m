function calibratedSpectra = calibrate_generic(rawSpectra,log,retrievalTool,TCold,calType)
%CALIBRATE_GENERIC Summary of this function goes here
%   Detailed explanation goes here

% Calibration type
%calType='time';

% CalibrationTime in Minute
calibrationTime=retrievalTool.calibrationTime;

% We define a complete cycle as 2-1-0-0-1-2:
hotIndCompleteCycle=find_completed_cycle(log,retrievalTool);
%return

% Group all indices by type for all the completed cycles:
ih=reshape([hotIndCompleteCycle; hotIndCompleteCycle+5],[],1);
ia=reshape([hotIndCompleteCycle+1;hotIndCompleteCycle+4],[],1);
ic=reshape([hotIndCompleteCycle+2;hotIndCompleteCycle+3],[],1);

% Convert it to int may help for the following
%hotIndCompleteCycle=uint32(hotIndCompleteCycle);

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
end

% Starting time for the complete cycle (will be the basis for separating
% the cycles according to calibrationTime). 
startingTimes=log.t(validStartIndices);

calibratedSpectra=struct();

switch calType
    case 'time'
        timeThresh=0:calibrationTime/60:24;
        
        % Storing the indices specific to each calibration cycle in a new
        % structure because by separating by time, we do not have the same
        % number of individual cycle per calibration cycle
        indices=struct();
        for i = 1:length(timeThresh)-1
            cond=startingTimes>timeThresh(i) & startingTimes<timeThresh(i+1);
            indice=[validStartIndices(cond); validStartIndices(cond)+1; validStartIndices(cond)+2; validStartIndices(cond)+3; validStartIndices(cond)+4; validStartIndices(cond)+5];
            indices(i).ind=indice;
        end
        cond=startingTimes>timeThresh(length(timeThresh)-1);
        lastIndices=[validStartIndices(cond); validStartIndices(cond)+1; validStartIndices(cond)+2; validStartIndices(cond)+3; validStartIndices(cond)+4; validStartIndices(cond)+5];
        indices(length(timeThresh)).ind=lastIndices;
        
        nCalibrationCycles=length(indices);
        % Based on the starting times, we will then group the cycles together to
        % form x minutes calibrationTime
        %
        rsAntenna=ones(nCalibrationCycles,retrievalTool.numberOfChannels)*NaN;
        rsCold=ones(nCalibrationCycles,retrievalTool.numberOfChannels)*NaN;
        rsHot=ones(nCalibrationCycles,retrievalTool.numberOfChannels)*NaN;
        THot=ones(nCalibrationCycles,1)*NaN;
        stdTHot=ones(nCalibrationCycles,1)*NaN;
        
        % We need to loop through the cycles because the number of averaged
        % spectra might be different between each cycles.
        for i=1:nCalibrationCycles
            ih=reshape(indices(i).ind([1,6],:),[],1);
            ia=reshape(indices(i).ind([2,5],:),[],1);
            ic=reshape(indices(i).ind([3,4],:),[],1);
            
            % Computing some useful quantities for this cycles:
            % Mean hot temperature for this cycle as well as its standard
            % deviation:
            THot(i)=nanmean(100*log.AI_0(reshape(indices(i).ind,[],1))); 
            stdTHot(i)=nanstd(100*log.AI_0(reshape(indices(i).ind,[],1)));
            
            % Number of hot/cold/antenna averaged spectra for this cycle
            % Considering all spectra that are not 100% NaN ... 
            nAvgSpectraHot(i)=length(ih)-sum(all(isnan(rawSpectra(ih,:)),1));
            nAvgSpectraAntenna(i)=length(ia)-sum(all(isnan(rawSpectra(ia,:)),1));
            nAvgSpectraCold(i)=length(ic)-sum(all(isnan(rawSpectra(ic,:)),1));
            
            % Mean and stdDev Temperature of the system for this cycle
            Tsys(i)=nanmean(log.FE_T_Sys(reshape(indices(i).ind,[],1)));
            stdTSys(i)=nanstd(log.FE_T_Sys(reshape(indices(i).ind,[],1)));
            
            % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
            rsHot(i,:)=nanmean(rawSpectra(ih,:),1);
            rsAntenna(i,:)=nanmean(rawSpectra(ia,:),1);
            rsCold(i,:)=nanmean(rawSpectra(ic,:),1);
        end
        
        % Calibration
        Tb = TCold + (THot-TCold).*(rsAntenna-rsCold)./(rsHot-rsCold);
        
    case 'all'
        % If no calibration time is provided, we calibrate every cycle
        % (2-1-0-0-1-2)
        indices=[validStartIndices; validStartIndices+1; validStartIndices+2; validStartIndices+3; validStartIndices+4; validStartIndices+5];
        % Number of calibration cycle for this day:
        nCalibrationCycles=size(indices,2);
        
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
        THot=nanmean(100*log.AI_0(indices),1);
        
        % Mean cold counts for each cycle (nCalibrationCycles x #channels)
        rsCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
        
        % Calibration
        Tb = TCold + (THot-TCold)'.*(rsAntenna-rsCold)./(rsHot-rsCold);
        
        % Mean System Temperature for each cycle
        Tsys=nanmean(log.FE_T_Sys(indices),1);
        % Std deviation of System Temperature for each cycle
        stdTSys=nanstd(log.FE_T_Sys(indices),1);
end

% Filling the structure as following:
% 
calibratedSpectra.Tb=Tb;
calibratedSpectra.nCalibrationCycles=nCalibrationCycles;
calibratedSpectra.Tsys=Tsys;

% Plotting N calibrated spectra to test
N=20;
l=floor(linspace(1,size(Tb,1),N));
% 
% figure();
% for i=1:N
%     plot(Tb(l(i),:))
%     ylim([100,370])
%     hold on
% end

    
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
%     function indices = find_indice(log,type,retrievalTool)
%         switch type
%             case 'cold'
%                 ind=retrievalTool.indiceCold;
%             case 'hot'
%                 ind=retrievalTool.indiceHot;
%             case 'antenna' 
%                 ind=retrievalTool.indiceAntenna;
%             otherwise
%             error('No valid type of indices provided')
%         end
%         indices=find(log.Position==ind & log.Tipping_Curve_active==0);
%     end   

end

