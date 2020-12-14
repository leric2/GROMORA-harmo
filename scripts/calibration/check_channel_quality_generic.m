function spectra = checking_channel_quality_gromos(spectra,calibrationTool,filterN)
%==========================================================================
% NAME          | checking_channel_quality_gromos.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 02.20
%               |
% ABSTRACT      | Check the channel quality of a spectrum (calibrated or
%               | integrated) depeding on the filter provided. 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:   1. spectra: calibrated or integratedSpectra to
%               |               check for channel quality. 
%               |           2. calibrationTool:
%               |           	- numberOfChannels
%               |           	- badChannels
%               |           	- filter1, filter2, ...
%               |           	- maxStdDevTbCal (optional)
%               |           	- maxStdDevTbInt (optional)
%               |           3. filterN: indice for the type of filtering to
%               |               perform
%               |
%               | OUTPUTS:  1. calibratedSpectra: completed with field channelsQuality
%               |
%==========================================================================
% Initially: all good
indicesGood=ones(calibrationTool.numberOfChannels,1);

% Bad channels for all calibrated cycles (dependant on the instruments):
indicesGood(calibrationTool.badChannels)=NaN;

% Defining the filering type to apply:
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
        % using the std dev of Tb to qualify good and bad spectra
        % (Calibrated spectra)
        for t = 1:length(spectra)
            indCyclesGood=indicesGood;
            
            indCyclesGood(spectra(t).stdTb > calibrationTool.maxStdDevTbCal) = NaN;
            
            spectra(t).channelsQuality=indCyclesGood';
        end
        return
    case 4
        % using the std dev of Tb to qualify good and bad spectra
        % (Integrated spectra)
        for t = 1:length(spectra)
            indCyclesGood=indicesGood;
            
            indCyclesGood(spectra(t).stdTb > calibrationTool.maxStdDevTbInt) = NaN;
            
            spectra(t).channelsQuality=indCyclesGood';
        end
        return
end

% creating boxcar filter (for the cases where we filter so).
boxCarFilter=ones(boxCarSize,1)/boxCarSize;

for t = 1:length(spectra)
    indCyclesGood=indicesGood;
    
    indCyclesGood(spectra(t).Tb>TbMax | spectra(t).Tb<TbMin)=NaN;
    
    TbFiltered=conv(spectra(t).Tb,boxCarFilter,'same');
    
    indCyclesGood(abs(spectra(t).Tb-TbFiltered)>boxCarthresh)=NaN;
    
    spectra(t).channelsQuality=indCyclesGood';
end
end