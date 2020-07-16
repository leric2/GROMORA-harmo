function calibratedSpectra = checking_channel_quality_gromos(calibratedSpectra,calibrationTool,filterN)
%==========================================================================
% NAME          | 
% TYPE          |
% AUTHOR(S)     |
% CREATION      |
%               |
% ABSTRACT      |
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
switch filterN
    case 1
        TbMax = calibrationTool.filter1.TbMax;
        TbMin= calibrationTool.filter1.TbMin;
        boxCarSize = calibrationTool.filter1.boxCarSize;
        boxCarthresh = calibrationTool.filter1.boxCarThresh;
    case 2
        TbMax = calibrationTool.filter2.TbMax;
        TbMin= calibrationTool.filter2.TbMin;
        boxCarSize = calibrationTool.filter2.boxCarSize;
        boxCarthresh = calibrationTool.filter2.boxCarThresh;
end

% creating boxcar filter
boxCarFilter=ones(boxCarSize,1)/boxCarSize;

indicesGood=ones(calibrationTool.numberOfChannels,1);

% Bad channels for all calibrated cycles (dependant on the instruments):
indicesGood(calibrationTool.badChannels)=NaN;

for t = 1:length(calibratedSpectra)
    indCyclesGood=indicesGood;
    
    indCyclesGood(calibratedSpectra(t).Tb>TbMax | calibratedSpectra(t).Tb<TbMin)=NaN;
    
    TbFiltered=conv(calibratedSpectra(t).Tb,boxCarFilter,'same');
    
    indCyclesGood(abs(calibratedSpectra(t).Tb-TbFiltered)>boxCarthresh)=NaN;
    
    calibratedSpectra(t).channelsQuality=indCyclesGood';
    
end


end