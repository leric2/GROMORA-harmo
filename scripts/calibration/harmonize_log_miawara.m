function logFile = harmonize_log_miawara(calibrationTool, logFile)
%==========================================================================
% NAME          | harmonize_log_miawara(logFile)
% TYPE          | function
% AUTHOR(S)     | Alistair Bell from modified code by F Schranz
% CREATION      | 2022-07-10
%               |
% ABSTRACT      | Modifications of the log file specific to MIAWARA
%               | 
%               |
% ARGUMENTS     | INPUTS: calibrationTool, logFile
%               |
%               | OUTPUTS: logFile (modified)
%               |
% CALLED by     | run_calibration(calibrationTool, logFile)
%               | import_MIAWARA_calibrationTool(calibrationTool)  
%               |
% CALLS         |
%               |
%               |

%==========================================================================
%
% Add variable time
logFile.time = datenum(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);

% Rainhood
logFile.rainhood_open = logFile.RH_open;
logFile.rainhood_closed = logFile.RH_closed;

%T sensors - Only 7 in hot load!
logFile.THot0=logFile.T_Hot_0;
logFile.THot1=logFile.T_Hot_1;
logFile.THot2=logFile.T_Hot_2;
logFile.THot3=logFile.T_Hot_3;
logFile.THot4=logFile.T_Hot_4;
logFile.THot5=logFile.T_Hot_5;
logFile.THot6=logFile.T_Hot_6;
logFile.THot7=logFile.T_Hot_7;

logFile.Mode_nr = logFile.meas_mode;
logFile.Mirror_pos = logFile.position;
disp('mirror pos:')
disp(all(~(logFile.position == 1)))

logFile.Elevation_Angle = (logFile.elevation_stop_deg_ + logFile.elevation_start_deg_ )/ 2;

%Calibration Modes
isTipping   = logFile.Mode_nr == 1;
isCustom    = logFile.Mode_nr == 2;
isLN2       = logFile.Mode_nr == 3;
isBalancing = logFile.Mode_nr == 6;

%Quality checks
rainhoodOK  = logFile.rainhood_open | logFile.rainhood_open == logFile.RH_closed;
angleOK     = ( logFile.Elevation_Angle < 85  &  logFile.Elevation_Angle > 5 );

% indices for calibration
logFile.isLine      = rainhoodOK & angleOK & (logFile.Mirror_pos == 2 | ...
    logFile.Mirror_pos == 5 | logFile.Mirror_pos == 7 ) &  (isBalancing  | isCustom) ;

logFile.isHot    = logFile.Mirror_pos == 0 & (isTipping  | isCustom | isLN2);
logFile.isColdSky = logFile.Mirror_pos == 4 & (isTipping  | isCustom | isLN2);
logFile.isRef     = rainhoodOK & logFile.Mirror_pos == 6 & (isBalancing | isCustom) ;
logFile.isColdLN2 = logFile.Mirror_pos == 1 & (isLN2 | isCustom);

%AB for investigating standing waves, the difference between the two mirror
%displacements are here made into a further category 
logFile.isMirrorDisp1 = logFile.x_stop_mm_ == 0;
logFile.isMirrorDisp2 = logFile.x_stop_mm_ == 3.3710;
%should probably add a check here to make sure that both these fields are
%filled - i.e. that the mirror displacement positions have not been changed

% name adaptations in case of polarisation and direction
logFile.dir = 0;
%logFile.dir = (logFile.Mode_nr == 11 |logFile.Mode_nr ==12) + 2*(logFile.Mode_nr==15 | logFile.Mode_nr == 16);

end
