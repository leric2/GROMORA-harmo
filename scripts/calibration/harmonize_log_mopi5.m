function logFile = harmonize_log_mopi5(calibrationTool, logFile)
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
% For MOPI5
if length(logFile.Year)==1
    logFile.Year=ones(length(logFile.t),1)*logFile.Year;
end
if length(logFile.Month)==1
    logFile.Month=ones(length(logFile.t),1)*logFile.Month;
end
if length(logFile.Day)==1
    logFile.Day=ones(length(logFile.t),1)*logFile.Day;
end

logFile.time = datenum(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);
logFile.dateTime = datetime(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);


logFile.Tipping_Curve_active=zeros(length(logFile.t),1);

% Hot temperature is in °C:
logFile.T_Hot_Absorber=logFile.T(:,1);

logFile.FE_T_Sys=-9999*ones(length(logFile.t),1);

logFile.FFT_adc_range=ones(length(logFile.t),1);
logFile.FFT_adc_overload=zeros(length(logFile.t),1);

logFile.Position=logFile.Mode_Target';
logFile.Elevation_Angle=logFile.Measurement_Elevation;

if ~(mean(logFile.Elevation_Angle(logFile.Position == 1)) == calibrationTool.elevationAngleHot)
    error('angle for the hot load might be wrongly defined') 
end

% log.T_Ceiling=log.TExt0;
% log.T_Floor=log.TExt1;
% log.T_Aircon_Out=log.TExt2;
% log.T_Window=log.TExt3;
% log.T_Amp1=log.TExt4;
% log.T_Amp2=log.TExt5;
% log.T_Mirror_View=log.TExt6;
% log.T_Reserved=log.TExt7;

logFile.LN2_Sensors_OK = ones(length(logFile.t),1);

logFile.LN2_Level_OK = ones(length(logFile.t),1);
    
logFile.T_Room=logFile.T(:,2);

end
