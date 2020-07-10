function log = harmonize_log_gromos(log)
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
% For GROMOS
% Add variable time
log.time = datenum(log.Year,log.Month,log.Day,log.Hour,log.Minute,log.Second);

% Hot temperature is in °C:
log.T_Hot_Absorber=log.T_Hot + 273.15; % to replace with calibrationTool.zeroDegInKelvin ?

log.T_Ceiling=log.TExt0;
log.T_Floor=log.TExt1;
log.T_Aircon_Out=log.TExt2;
log.T_Window=log.TExt3;
log.T_Amp1=log.TExt4;
log.T_Amp2=log.TExt5;
log.T_Mirror_View=log.TExt6;
log.T_Reserved=log.TExt7;

log.LN2_Sensors_OK = ~log.LN2_Relay;

log.LN2_Level_OK = log.LN2_above_High & ~(log.LN2_above_Low);
    
% log.T_Room;

end
