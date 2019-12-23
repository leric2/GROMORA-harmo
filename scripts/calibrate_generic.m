function calibratedSpectra = calibrate_generic(rawSpectra,log,retrievalTool,TCold)
%CALIBRATE_GENERIC Summary of this function goes here
%   Detailed explanation goes here

% For the beginning we just group a number of given cycle together

% CalibrationTime in Minute
calibrationTime=10;

% We define a complete cycle as 2-1-0-0-1-2:
hotIndCompleteCycle=find_completed_cycle(log,retrievalTool);

% Convert it to int may help for the following
hotIndCompleteCycle=uint32(hotIndCompleteCycle);

% Starting time for the complete cycle
startingTimes=log.t(hotIndCompleteCycle);

timeThresh=0:calibrationTime/60:24;

data=struct();
for i = 1:length(timeThresh)-1
    cond=startingTimes>timeThresh(i) & startingTimes<timeThresh(i+1);
    indices=[hotIndCompleteCycle(cond); hotIndCompleteCycle(cond)+1; hotIndCompleteCycle(cond)+2; hotIndCompleteCycle(cond)+3; hotIndCompleteCycle(cond)+4; hotIndCompleteCycle(cond)+5];
    %data.(strcat('Calib',num2str(i)))=indices;
    data(i).ind=indices;
end
cond=startingTimes>timeThresh(length(timeThresh)-1);
lastIndices=[hotIndCompleteCycle(cond); hotIndCompleteCycle(cond)+1; hotIndCompleteCycle(cond)+2; hotIndCompleteCycle(cond)+3; hotIndCompleteCycle(cond)+4; hotIndCompleteCycle(cond)+5];
data(length(timeThresh)).ind=lastIndices;

% Based on the starting times, we will then group the cycles together to
% form x minutes calibrationTime
%
% indices will then be a 3D matrix ?? (#FinalCycle x 6 x #ofAppendedCycles)
%

% If no calibration time is provided, we calibrate every cycle
% (2-1-0-0-1-2)
if ~isfield(retrievalTool,'calibrationTime')
    indices=[hotIndCompleteCycle; hotIndCompleteCycle+1; hotIndCompleteCycle+2; hotIndCompleteCycle+3; hotIndCompleteCycle+4; hotIndCompleteCycle+5];
    ih=[hotIndCompleteCycle;hotIndCompleteCycle+5];
    ia=[hotIndCompleteCycle+1;hotIndCompleteCycle+4];
    ic=[hotIndCompleteCycle+2;hotIndCompleteCycle+3];
    
    % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
    rSAntenna=nanmean(cat(3,rawSpectra(ia(1,:),:),rawSpectra(ia(2,:),:)),3);
    
    % Mean hot counts for each cycle (nCalibrationCycles x #channels)
    rSHot=nanmean(cat(3,rawSpectra(ih(1,:),:),rawSpectra(ih(2,:),:)),3);
    THot=nanmean(100*log.AI_0(indices),1);
    
    % Mean cold counts for each cycle (nCalibrationCycles x #channels)
    rSCold=nanmean(cat(3,rawSpectra(ic(1,:),:),rawSpectra(ic(2,:),:)),3);
else
    for i=1:length(data)
        ih=data(i).ind([1,6],:);
        ia=data(i).ind([2,5],:);
        ic=data(i).ind([3,4],:);
        % Mean Antenna counts for each cycle (nCalibrationCycles x #channels)
        
    end

end
% Number of calibration cycle for this day:
nCalibrationCycles=size(indices,2);

% We could loop into the different calibration cycles but we will try to
% avoid it and perform all the operation for all calibration cycles
% simultaneously
%rSpect=ones(6,nCalibrationCycles,retrievalTool.numberOfChannels);




% Calibration
Tb=TCold + (THot-TCold)'.*(rSAntenna-rSCold)./(rSHot-rSCold);




%     % Raw spectra corresponding to each completed cycle
%     rSHot1=rawSpectra(hotIndCompleteCycle,:);   rSHot2=rawSpectra(hotIndCompleteCycle+5,:);
%     rSCold1=rawSpectra(hotIndCompleteCycle+2,:); rSCold2=rawSpectra(hotIndCompleteCycle+3,:);
%     rSAntenna1=rawSpectra(hotIndCompleteCycle+1,:); rSAntenna2=rawSpectra(hotIndCompleteCycle+4,:);
%
%     % plot(rawSpectra([5 ],:))
%
%
%     rSHot=mean(cat(3,rSHot1,rSHot2),3);
%     rSCold=mean(cat(3,rSCold1,rSCold2),3);
%     rSAntenna=mean(cat(3,rSAntenna1,rSAntenna2),3);


% Mean System Temperature for each cycle
Tsys=nanmean(log.FE_T_Sys(indices),1);

% Std deviation of System Temperature for each cycle
stdTSys=nanstd(log.FE_T_Sys(indices),1);


% The actual calibration
% Some metrics



    
    function firstIndCompleteCycle = find_completed_cycle(log,retrievalTool)
        indCold=retrievalTool.indiceCold;
        indHot=retrievalTool.indiceHot;
        indAntenna=retrievalTool.indiceAntenna;
        hotInd=find(log.Tipping_Curve_active==0 & (log.Position==indHot));
        extendedPos=[log.Position -9999 -9999 -9999 -9999 -9999];
        firstIndCompleteCycle=hotInd((extendedPos(hotInd+1)==indAntenna & extendedPos(hotInd+2)==indCold & extendedPos(hotInd+3)==indCold & extendedPos(hotInd+4)==indAntenna & extendedPos(hotInd+5)==indHot));
    end
    
    
    
    
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

