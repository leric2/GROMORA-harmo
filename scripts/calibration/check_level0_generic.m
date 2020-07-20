function warningLevel0 = check_level0_generic(logFile,rawSpectra,calibrationTool)
%==========================================================================
% NAME          | check_level0_generic 
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
% Initialize structure containing the error that are non fatal for the
% retrieval
warningLevel0='';

% Compare size of the log and raw spectra for this day
mLog = size(logFile.x,2);

% Check if the size of the binary is consistent with the log
w=whos('rawSpectra');
try
    assert(w.bytes()==(mLog*calibrationTool.numberOfChannels*calibrationTool.bytesPerValue),'consistency:unconsistentBinarySize_','Size of the binary unconsistent');
catch ME
    warningLevel0=append(warningLevel0,ME.identifier);
end

% Check the number of channels
% try
%     assert(n==retrievalTool.numberOfChannels,'consistency:numberOfChannel');
% catch ME
%     errorLevel0_1a.nChannels=ME.identifier;
% end

% Number of tipping curve calibration for this day
nTipping=sum(logFile.Tipping_Curve_active==1)/calibrationTool.tippingSize;

try
    assert(((nTipping<calibrationTool.numberOfTippingCurveExpected+calibrationTool.toleranceTippingCurves)&&(nTipping>calibrationTool.numberOfTippingCurveExpected-calibrationTool.toleranceTippingCurves)) ,'consistency:numberTippingCurve_','Number of Tipping unconsistent');
catch ME
    % warningLevel0=append(warningLevel0,ME.identifier);
end

% Approximate number of cycles completed for this day
nCycles=sum(logFile.Tipping_Curve_active==0)/6;
try
    assert(((nCycles<calibrationTool.numberOfCyclesExpected+calibrationTool.toleranceNumberCycles)&&(nCycles>calibrationTool.numberOfCyclesExpected-calibrationTool.toleranceNumberCycles)),'consistency:numberCycles_','Number of cycles unconsistent');
catch ME
    warningLevel0=append(warningLevel0,ME.identifier);
end

% Temporal check to keep only exact daily file for the calibration
t0 = datenum(calibrationTool.dateStr,'yyyy_mm_dd');
isDayAfter = logFile.time >= t0+1;
isDayBefore = logFile.time < t0;
isExtra = (isDayAfter | isDayBefore);

if sum(isExtra) > 0    
    f = [calibrationTool.extraFileFolder calibrationTool.filename '_extra'];
    
    M = rawSpectra(isExtra,:)';
    fid = fopen([f '.bin'], 'w');
    fwrite(fid,M(:));
    fclose(fid);
    
    writecell(logFile.header',[f '.txt'],'Delimiter',calibrationTool.delimiter_logfile)
    dlmwrite([f '.txt'],logFile.x(:,isExtra)','delimiter',calibrationTool.delimiter_logfile,'-append');
    disp('Some extra timestamp have been saved, please make sure to use "read_level0_missing" function');
    warningLevel0=append(warningLevel0,'extra_timestamp');
end
end

