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
indicesGood=ones(calibrationTool.numberOfChannels,1);

% Bad channels for all calibrated cycles (dependant on the instruments):
indicesGood(calibrationTool.badChannels)=NaN;

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
    case 3
        % using the std dev of Tb to quality good and bad spectra
        for t = 1:length(calibratedSpectra)
            indCyclesGood=indicesGood;
            
            indCyclesGood(calibratedSpectra(t).stdTb > calibrationTool.maxStdDevTbCal) = NaN;
            
            calibratedSpectra(t).channelsQuality=indCyclesGood';
        end
        return
    case 4
        % using the std dev of Tb to quality good and bad spectra (for
        % Integrated spectra
        for t = 1:length(calibratedSpectra)
            indCyclesGood=indicesGood;
            
            indCyclesGood(calibratedSpectra(t).stdTb > calibrationTool.maxStdDevTbInt) = NaN;
            
            calibratedSpectra(t).channelsQuality=indCyclesGood';
        end
        return
end

% creating boxcar filter
boxCarFilter=ones(boxCarSize,1)/boxCarSize;

for t = 1:length(calibratedSpectra)
    indCyclesGood=indicesGood;
    
    indCyclesGood(calibratedSpectra(t).Tb>TbMax | calibratedSpectra(t).Tb<TbMin)=NaN;
    
    TbFiltered=conv(calibratedSpectra(t).Tb,boxCarFilter,'same');
    
    indCyclesGood(abs(calibratedSpectra(t).Tb-TbFiltered)>boxCarthresh)=NaN;
    
    calibratedSpectra(t).channelsQuality=indCyclesGood';
    
end


end