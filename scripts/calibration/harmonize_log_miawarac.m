

function log = harmonize_log_miawarac(log)
%==========================================================================
% NAME          | harmonize_log_miawarac(log)
% TYPE          | function
% AUTHOR(S)     | Franziska Schranz
% CREATION      | 2020-06-08
%               |
% ABSTRACT      | Modifications of the log file specific to MIAWARA-C
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: log
%               |
%               | OUTPUTS: log (modified)
%               |
% CALLED by     | run_calibration(calibrationTool)
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
% For MIAWARA-C

% Add variable time
log.time = datenum(log.Year,log.Month,log.Day,log.Hour,log.Minute,log.Second)

% Rainhood
g = log.DIO1;
openbit= 4; 
closedbit = 2;
log.rainhood_open = (bitand(g,openbit))==0;
log.rainhood_closed = (bitand(g,closedbit))==0;

% Nyalesund: Bad contact, TH5-Th7 = copy of THot0-4
if datenum(log.time(1,:),31)>datenum(2015,09,01)
log.THot4=log.THot0;
log.THot5=log.THot1;
log.THot6=log.THot2;
log.THot7=log.THot3;
end

% after June 2013: THot8-15 = copy of THot0-7
if datenum(log.time(1,:),31)>datenum(2013,06,01)
log.THot8=log.THot0;
log.THot9=log.THot1;
log.THot10=log.THot2;
log.THot11=log.THot3;
log.THot12=log.THot4;
log.THot13=log.THot5;
log.THot14=log.THot6;
log.THot15=log.THot7;
end




end
