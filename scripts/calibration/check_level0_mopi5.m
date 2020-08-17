function warningLevel0 = check_level0_mopi5(logFile,rawSpectra,calibrationTool)
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

% Temporal check to keep only exact daily file for the calibration
t0 = datenum(calibrationTool.dateStr,'yyyy_mm_dd');
isDayAfter = logFile.time >= t0+1;
isDayBefore = logFile.time < t0;
isExtra = (isDayAfter | isDayBefore);

if sum(isExtra) > 0    
    f = [calibrationTool.extraFileFolder calibrationTool.filename '_extra'];
    
    M = rawSpectra(isExtra,:)';
    fid = fopen([f '_' calibrationTool.spectrometer calibrationTool.binaryDataExtension], 'w');
    fwrite(fid,M(:));
    fclose(fid);
    
    writecell(logFile.header,[f calibrationTool.logFileDataExtension],'Delimiter',calibrationTool.delimiter_logfile)
    dlmwrite([f calibrationTool.logFileDataExtension], logFile.x(isExtra,:),'delimiter',calibrationTool.delimiter_logfile,'-append');
    disp('Some extra timestamp have been saved, please make sure to use "read_level0_missing" function');
    warningLevel0=append(warningLevel0,'extra_timestamp');
end
end

