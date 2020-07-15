function logFile = harmonize_log_miawarac(calibrationTool, logFile)
%==========================================================================
% NAME          | harmonize_log_miawarac(logFile)
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
logFile.time = datenum(logFile.Year,logFile.Month,logFile.Day,logFile.Hour,logFile.Minute,logFile.Second);

% Rainhood
g = logFile.DIO1;
openbit= 4; 
closedbit = 2;
logFile.rainhood_open = (bitand(g,openbit))==0;
logFile.rainhood_closed = (bitand(g,closedbit))==0;

% Nyalesund: Bad contact, TH5-Th7 = copy of THot0-4
if logFile.time(1,:)>datenum(2015,09,01)
logFile.THot4=logFile.THot0;
logFile.THot5=logFile.THot1;
logFile.THot6=logFile.THot2;
logFile.THot7=logFile.THot3;
end

% after June 2013: THot8-15 = copy of THot0-7
if logFile.time(1,:)>datenum(2013,06,01)
logFile.THot8=logFile.THot0;
logFile.THot9=logFile.THot1;
logFile.THot10=logFile.THot2;
logFile.THot11=logFile.THot3;
logFile.THot12=logFile.THot4;
logFile.THot13=logFile.THot5;
logFile.THot14=logFile.THot6;
logFile.THot15=logFile.THot7;
end

% % add additional field with status(k) balancing
% for k = 1:length(logFile.Year)
%     if logFile.Mode_nr(k) == 0
%         status(k)= 2;   %'tipping_all_sky'; 
%     elseif logFile.Mode_nr(k) == 2
%         fprintf('ALLAN VARIANCE => not in db, %s \n',logFile.time(k,:)); 
%         k=k+1;
%         continue;
%     elseif logFile.Mode_nr(k) == 3
%         status(k)= 3;   %'tipping';
%     elseif logFile.Mode_nr(k) == 8
%         status(k)= 4;   %'calib_ln2';
%     elseif logFile.Mode_nr(k) == 9
%         status(k)= 5;   %'sun_scan';
%     elseif logFile.Mode_nr(k) == 10
%         status(k)= 1;   %'balancing'; %search reference
%         balancing_mode=3;
%     elseif logFile.Mode_nr(k) == 11 && logFile.program(k) == 1
%         status(k)= 1;   %'balancing'; %search sky, balancing left
%         balancing_mode=1;
%     elseif logFile.Mode_nr(k) == 12 && logFile.program(k) == 2
%         status(k)= 1;   %'balancing'; %search sky with auto-ref, balancing left
%         balancing_mode=2;
%     elseif logFile.Mode_nr(k) == 15 && logFile.program(k) == 5
%         status(k)= 1;   %'balancing'; %search sky, balancing right
%         balancing_mode=1;
%     elseif logFile.Mode_nr(k) == 16 && logFile.program(k) == 6
%         status(k)= 1;   %'balancing'; %search sky with auto-ref, balancing right
%         balancing_mode=2;  
% 
%         %%% dual balancing
%     elseif logFile.Mode_nr(k) == 11 && logFile.program(k) == 7
%         status(k)= 1;   %'dual_balancing_left'; %dual balancing, search sky, balancing left
%         balancing_mode=1;
%     elseif logFile.Mode_nr(k) == 12 && logFile.program(k) == 8
%         status(k)= 1;   %'dual_balancing_left'; %dual balancing, search sky with auto-ref, balancing left
%         balancing_mode=2;
%     elseif logFile.Mode_nr(k) == 15 && logFile.program(k) == 7
%         status(k)= 1;   %'dual_balancing_right'; %dual balancing, search sky, balancing right
%         balancing_mode=1;
%     elseif logFile.Mode_nr(k) == 16 && logFile.program(k) == 8
%         status(k)= 1;   %'dual_balancing_right'; %dual balancing, search sky with auto-ref, balancing right
%         balancing_mode=2;  
%     end
% end
  


% indices for calibration

isBalancing = logFile.Mode_nr == 11 | logFile.Mode_nr == 12 | logFile.Mode_nr == 15 | logFile.Mode_nr == 16 ;
isTipping   = logFile.Mode_nr == 3;
rainhoodOK  = logFile.rainhood_open | logFile.rainhood_open == logFile.RH_closed;
angleOK     = ( logFile.Mirror_elevation < 65 & logFile.Mirror_elevation  >  5 ) | ( logFile.Mirror_elevation < 175 & logFile.Mirror_elevation > 140 ); % for line


logFile.isLine      = isBalancing & rainhoodOK & angleOK & (logFile.Mirror_pos == 0 | logFile.Mirror_pos == 5 | logFile.Mirror_pos == 6 | logFile.Mirror_pos == 7);
logFile.isHot       = logFile.Mirror_pos == 1 & isTipping;
logFile.isColdSky   = logFile.Mirror_pos == 3 & isTipping;
logFile.isRef       = isBalancing & rainhoodOK & logFile.Mirror_pos == 4;


% name adaptations

logFile.Elevation_Angle = logFile.Mirror_elevation;

logFile.dir = (logFile.Mode_nr == 11 |logFile.Mode_nr ==12) + 2*(logFile.Mode_nr==15 | logFile.Mode_nr == 16);

% %%% dual balancing
%         elseif log1.Mode_nr(k) == 11 && log1.program(k) == 7
%             status='dual_balancing_left'; %dual balancing, search sky, balancing left
%             balancing_mode=1;
%         elseif log1.Mode_nr(k) == 12 && log1.program(k) == 8
%             status='dual_balancing_left'; %dual balancing, search sky with auto-ref, balancing left
%             balancing_mode=2;
%         elseif log1.Mode_nr(k) == 15 && log1.program(k) == 7
%             status='dual_balancing_right'; %dual balancing, search sky, balancing right
%             balancing_mode=1;
%         elseif log1.Mode_nr(k) == 16 && log1.program(k) == 8
%             status='dual_balancing_right'; %dual balancing, search sky with auto-ref, balancing right
%             balancing_mode=2;  
%     end


end
