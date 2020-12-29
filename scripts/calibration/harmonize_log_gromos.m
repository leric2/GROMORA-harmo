function logFile = harmonize_log_gromos(calibrationTool,logFile)
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
%   Spectr_left_wing_start	Spectr_left_wing_width	Spectr_line_center	Spectr_line_width	Spectr_T_Line_Amp	Spectr_T_Peak	Spectr_T_Wing	
%   Data_file_size	SW_version	IWV
%   
%==========================================================================

% For GROMOS
% Add variable time
logFile.time = datenum(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);
logFile.dateTime = datetime(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second, 'TimeZone',calibrationTool.timeZone);

if (~(mean(logFile.Elevation_Angle(logFile.Position == 0)) == calibrationTool.elevationAngleCold) && calibrationTool.dateTime>datetime(2012,04,26, 'TimeZone',calibrationTool.timeZone))
    error('angle for the cold load might be wrongly defined') 
end

% Hot temperature is in °C:
logFile.T_Hot_Absorber=logFile.T_Hot + calibrationTool.zeroDegInKelvin; % to replace with  ?

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
    logFile.T_Ceiling=logFile.AI_0;
    logFile.T_Floor=logFile.AI_0;
    logFile.T_Aircon_Out=logFile.AI_0;
    logFile.T_Window=logFile.AI_0;
    logFile.T_Amp1=logFile.AI_0;
    logFile.T_Amp2=logFile.AI_0;
    logFile.T_Mirror_View=logFile.AI_0;
    logFile.T_Reserved=logFile.AI_0;
end
logFile.LN2_Sensors_OK = ~logFile.LN2_Relay;

logFile.LN2_Level_OK = (logFile.LN2_above_High == calibrationTool.goodFlagLN2Above) & (logFile.LN2_above_Low == calibrationTool.goodFlagLN2Below);
    
% log.T_Room;

end
