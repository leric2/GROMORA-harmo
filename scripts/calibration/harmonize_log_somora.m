function log = harmonize_log_somora(log)
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
% SOMORA
% Add variable time
log.time = datenum(log.Year,log.Month,log.Day,log.Hour,log.Minute,log.Second);
log.dateTime = datetime(log.Year,log.Month,log.Day,log.Hour,log.Minute,log.Second);

log.Tipping_Curve_active=(log.Position==6);
log.T_Hot_Absorber=log.AI_0*100;
log.T_Window=log.AI_1*100;
log.T_Out=log.AI_7*100;
log.T_Room=log.AI_1*100;
end
