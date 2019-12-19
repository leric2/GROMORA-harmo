function errorLevel0_1a = check_level0_generic(log,rawSpectra,retrievalTool,errorLevel0_1a)
% 
% Quality check for the raw spectra (level0 data)
% 
%
% INPUT
% log:          Housekeeping structure with field names for each parameter
% rawSpectra:   [N,channels] array of binary file
%
% OUTPUT
%
%
% DATA FORMAT
%
%
% Revisions

% Compare size of the log and raw spectra for this day
[m,n] = size(rawSpectra);
mLog = size(log.x,2);
try
    assert(mLog==m,'consistency:rawSpectraNotCorrespondingLog');
catch ME
    errorLevel0_1a.size=ME.identifier;
end

% Check the number of channels
try
    assert(n==retrievalTool.numberOfChannels,'consistency:numberOfChannel');
catch ME
    errorLevel0_1a.nChannels=ME.identifier;
end

% Number of tipping curve calibration for this day
nTipping=sum(log.Tipping_Curve_active==1)/retrievalTool.tippingSize;

try
    assert(((nTipping<retrievalTool.numberOfTippingCurveExpected+retrievalTool.toleranceTippingCurves)&&(nTipping>retrievalTool.numberOfTippingCurveExpected-retrievalTool.toleranceTippingCurves)) ,'consistency:numberTippingCurve');
catch ME
    errorLevel0_1a.nTipping=ME.identifier;
end

% Approximate number of cycles completed for this day
nCycles=sum(log.Tipping_Curve_active==0)/6;

try
    assert((nCycles<retrievalTool.numberOfCyclesExpected+retrievalTool.toleranceNumberCycles)&&(nCycles>retrievalTool.numberOfCyclesExpected-retrievalTool.toleranceNumberCycles),'consistency:numberCycles');
catch ME
    errorLevel0_1a.nCycles=ME.identifier;
end
end

