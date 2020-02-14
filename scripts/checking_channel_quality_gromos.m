function calibratedSpectra = checking_channel_quality_gromos(calibratedSpectra,retrievalTool)
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
retrievalTool.TbMax=300;

retrievalTool.TbMin=20;

retrievalTool.boxCarSize=51;
% retrievalTool.badChannels=[16384 16385];
retrievalTool.badChannels=[];
retrievalTool.boxCarThreshold=7;

% creating boxcar filter
boxCarFilter=ones(retrievalTool.boxCarSize,1)/retrievalTool.boxCarSize;


indicesGood=ones(retrievalTool.numberOfChannels,1);

% Bad channels for all calibrated cycles (dependant on the instruments):
indicesGood(retrievalTool.badChannels)=0;

for t = 1:length(calibratedSpectra)
    indCyclesGood=indicesGood;
    
    indCyclesGood(calibratedSpectra(t).Tb>retrievalTool.TbMax | calibratedSpectra(t).Tb<retrievalTool.TbMin)=0;
    
    TbFiltered=conv(calibratedSpectra(t).Tb,boxCarFilter,'same');
    
    indCyclesGood(abs(calibratedSpectra(t).Tb-TbFiltered)>retrievalTool.boxCarThreshold)=0;
    
    calibratedSpectra(t).channelsQuality=indCyclesGood;
    
end


end