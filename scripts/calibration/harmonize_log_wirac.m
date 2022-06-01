function logFile = harmonize_log_wirac(calibrationTool,logFile)
%==========================================================================
% NAME          | harmonize_log_gromos(calibrationTool,log)
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | harmonizing log of GROMOS to get standard log structure
%               | for the calibration.
%               |
% ARGUMENTS     | INPUTS: 1. calibrationTool:
%               |           - timeZone
%               |           - elevationAngleCold
%               |           - goodFlagLN2Below
%               |           - goodFlagLN2Above
%               |           - 
%               |           
%               |         2. logFile: original raw log file read
%               |
%               | OUTPUTS: - logFile: harmonized log file for this
%               |           instrument.
%               | 
%               |
%
% As standardized output variables, we want (all temperature in Kelvin):
%   T_Room
%   FE_T_Sys
%   Year Month Day Hour Minute Second
%   T_Hot_Absorber measured on the absorber
%   T_Window
%   T_Out
%   FE_T_Sys
%   Position Elevation_Angle	
%   Tipping_Curve_active Tipping_Angle_Nr	
%   Ferranti_Lock	PLL_Lock	V_Gunn	
%   LN2_above_High	LN2_above_Low	LN2_Relay	
%   FFT_adc_range	FFT_adc_overload	FFT_T_FPGA	FFT_Mode	FFT_Nr_of_acq	
%   
%==========================================================================

% For GROMOS
% Add variable time
logFile.time = datenum(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);
logFile.dateTime = datetime(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second, 'TimeZone',calibrationTool.timeZone);


% Hot temperature:
logFile.T_Hot_Absorber = (logFile.HousingParameters_Temp_6 + logFile.HousingParameters_Temp_7)/2

if isfield(logFile,'TExt0')
    logFile.T_Ceiling=logFile.TExt0;
    logFile.T_Floor=logFile.TExt1;
    logFile.T_Aircon_Out=logFile.TExt2;
    logFile.T_Window=logFile.TExt3;
    logFile.T_Amp1=logFile.TExt4;
    logFile.T_Amp2=logFile.TExt5;
    logFile.T_Mirror_View=logFile.TExt6;
    logFile.T_Reserved=logFile.TExt7;
else
    %TODO
    %disp('Error with log file parameters')
    % TOCHECK if this is really some kind of room temperature
   % logFile.T_Ceiling=logFile.AI_7*100;
%     logFile.T_Floor=logFile.AI_0;
%     logFile.T_Aircon_Out=logFile.AI_0;
%     logFile.T_Window=logFile.AI_0;
%     logFile.T_Amp1=logFile.AI_0;
%     logFile.T_Amp2=logFile.AI_0;
%     logFile.T_Mirror_View=logFile.AI_0;
%     logFile.T_Reserved=logFile.AI_0;
end
logFile.LN2_Sensors_OK = ones(length(logFile.t),1);

logFile.LN2_Level_OK =ones(length(logFile.t),1)

logFile.Freq_Lock = ones(length(logFile.t),1)

logFile.T_Room = nan*ones(length(logFile.t),1);
logFile.comment= "wirac measurement campaing"
if ~isfield(logFile, 'FE_T_Sys')
    logFile.FE_T_Sys = ones(length(logFile.t),1)*nan;
    
end
logFile.Position = logFile.Mode_Target_0;
logFile.Elevation_Angle = logFile.Measurement_Elevation
logFile.FFT_adc_overload = zeros(length(logFile.t),1)
logFile.Tipping_Curve_active = zeros(length(logFile.t),1);
if calibrationTool.timeNumber < datenum(2010,05,13)
    % Seems that before this date, there was no room temperature measured.
    logFile.T_Room = ones(length(logFile.t),1)*nan;
elseif (calibrationTool.timeNumber > datenum(2010,08,18) && calibrationTool.timeNumber < datenum(2012,07,23))
    %  Removing dependance on the Above_Low flags during this period.
    logFile.LN2_Level_OK = (logFile.LN2_above_High == calibrationTool.goodFlagLN2Above);
end
