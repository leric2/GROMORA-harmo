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

% Hot temperature is in °C:
log.T_Hot_Absorber=log.T_Hot+273.15;

% log.T_Room;

end
