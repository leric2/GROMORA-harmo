function logFile = harmonize_log_somora(calibrationTool, logFile)
%==========================================================================
% NAME          | harmonize_log_somora(calibrationTool,log)
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | harmonizing log of SOMORA to get standard log structure
%               | for the calibration.
%               |
%               |
% ARGUMENTS     | INPUTS: 1. calibrationTool:
%               |           - timeZone
%               |           - elevationAngleCold

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
%   
% SOMORA
% Add variable time
logFile.time = datenum(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);
logFile.dateTime = datetime(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second, 'TimeZone',calibrationTool.timeZone);

% SOMORA has no TC flags but an Position indice for it.
logFile.Tipping_Curve_active=(logFile.Position==calibrationTool.indiceTC);

% Some temperature sensors:
logFile.T_Hot_Absorber=logFile.AI_0*100;
logFile.T_Window=logFile.AI_3*100;
logFile.T_Out=logFile.AI_7*100;
logFile.T_Room=logFile.AI_1*100;

logFile.Freq_Lock = (logFile.IF_LO1_Lock & logFile.IF_LO2_Lock);   
% Sometimes the log file flags are single variable for the whole day
% so we convert it to vector
if length(logFile.Freq_Lock)==1
    logFile.Freq_Lock = ones(length(logFile.t),1)*logFile.Freq_Lock;
end
if length(logFile.FFT_Nr_of_acq)==1
    logFile.FFT_Nr_of_acq = ones(length(logFile.t),1)*logFile.FFT_Nr_of_acq;

end
if length(logFile.LN2_Level_OK)==1
    logFile.LN2_Level_OK = ones(length(logFile.t),1)*logFile.LN2_Level_OK;
end
if ~isfield(logFile, 'comment')
    logFile.comment = 'warning: corrupted raw file';
end
