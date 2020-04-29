function warningLevel0 = check_level0_generic(log,rawSpectra,retrievalTool)
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
% 
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

% Initialize structure containing the error that are non fatal for the
% retrieval
warningLevel0='';

% Compare size of the log and raw spectra for this day
mLog = size(log.x,2);

% Check if the size of the binary is consistent with the log
w=whos('rawSpectra');
try
    assert(w.bytes()==(mLog*retrievalTool.numberOfChannels*retrievalTool.bytesPerValue),'consistency:unconsistentBinarySize_','Size of the binary unconsistent');
catch ME
    warningLevel0=append(warningLevel0,ME.identifier);
end

% Check the number of channels
% n=
% try
%     assert(n==retrievalTool.numberOfChannels,'consistency:numberOfChannel');
% catch ME
%     errorLevel0_1a.nChannels=ME.identifier;
% end

% Number of tipping curve calibration for this day
nTipping=sum(log.Tipping_Curve_active==1)/retrievalTool.tippingSize;

try
    assert(((nTipping<retrievalTool.numberOfTippingCurveExpected+retrievalTool.toleranceTippingCurves)&&(nTipping>retrievalTool.numberOfTippingCurveExpected-retrievalTool.toleranceTippingCurves)) ,'consistency:numberTippingCurve_','Number of Tipping unconsistent');
catch ME
    warningLevel0=append(warningLevel0,ME.identifier);
end

% Approximate number of cycles completed for this day
nCycles=sum(log.Tipping_Curve_active==0)/6;
try
    assert(((nCycles<retrievalTool.numberOfCyclesExpected+retrievalTool.toleranceNumberCycles)&&(nCycles>retrievalTool.numberOfCyclesExpected-retrievalTool.toleranceNumberCycles)),'consistency:numberCycles_','Number of cycles unconsistent');
catch ME
    warningLevel0=append(warningLevel0,ME.identifier);
end

end

