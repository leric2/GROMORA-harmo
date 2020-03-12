function calibratedSpectra = checking_channel_quality_gromos(calibratedSpectra,retrievalTool,filterN)
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
        TbMax = retrievalTool.filter1.TbMax;
        TbMin= retrievalTool.filter1.TbMin;
        boxCarSize = retrievalTool.filter1.boxCarSize;
        boxCarthresh = retrievalTool.filter1.boxCarThresh;
    case 2
        TbMax = retrievalTool.filter2.TbMax;
        TbMin= retrievalTool.filter2.TbMin;
        boxCarSize = retrievalTool.filter2.boxCarSize;
        boxCarthresh = retrievalTool.filter2.boxCarThresh;
end

% creating boxcar filter
boxCarFilter=ones(boxCarSize,1)/boxCarSize;

indicesGood=ones(retrievalTool.numberOfChannels,1);

% Bad channels for all calibrated cycles (dependant on the instruments):
indicesGood(retrievalTool.badChannels)=NaN;

for t = 1:length(calibratedSpectra)
    indCyclesGood=indicesGood;
    
    indCyclesGood(calibratedSpectra(t).Tb>TbMax | calibratedSpectra(t).Tb<TbMin)=NaN;
    
    TbFiltered=conv(calibratedSpectra(t).Tb,boxCarFilter,'same');
    
    indCyclesGood(abs(calibratedSpectra(t).Tb-TbFiltered)>boxCarthresh)=NaN;
    
    calibratedSpectra(t).channelsQuality=indCyclesGood';
    
end


end