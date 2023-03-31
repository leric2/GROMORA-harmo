function [difference_hot, difference_reference] = differences_signal_positions(rawSpectra, calibrationTool, logFile)

%==========================================================================
% NAME          | differences_signal_positions.m
% TYPE          | function
% AUTHOR(S)     | Alistair Bell 
% CREATION      | 11.2022
%               |
% ABSTRACT      | Reading the raw input at different postitions to be
%               | see if standing waves are indeed a problem for the 
%               | radiometer
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: 1. rawSpectra: Matrix of raw data (1 line per
%               |           cycle, channels as columns.
%               |         2. logFile: standardized log file
%               |         3. calibrationTool:
%               |               - numberOfChannels
%               |               - TC
%               |               - indiceHot, indiceAntenna, indiceCold,
%               |                 indiceTC
%               |               - referenceTime
%               |               - tippingSize
%               |              
%               |
%               | OUTPUT: [pos1_signal, pos2_singal]
%               |        raw signal by the positions indicated
%               |
%==========================================================================
disp('Starting signal position Difference')
%Get indices of hot load observations at both positions

pos1_hot_ind = find(logFile.Mirror_pos == 0 & logFile.x_stop_mm_<1);
pos2_hot_ind = find(logFile.Mirror_pos == 0 & logFile.x_start_mm_>3);

s_hot_pos1 = rawSpectra (pos1_hot_ind,:);
s_hot_pos2 = rawSpectra (pos2_hot_ind,:);

difference_hot = (mean(s_hot_pos1, 1, "omitnan") - mean(s_hot_pos2, 1, "omitnan"))./mean(s_hot_pos1, 1, "omitnan");

%Get indices of reference load observations at both positions
pos1_ref_ind = find(logFile.Mirror_pos == 6 & logFile.x_stop_mm_<1);
pos2_ref_ind = find(logFile.Mirror_pos == 6 & logFile.x_start_mm_>3);

s_ref_pos1 = rawSpectra (pos1_ref_ind,:);
s_ref_pos2 = rawSpectra (pos2_ref_ind,:);

difference_reference = (mean(s_ref_pos1, 1, "omitnan") - mean(s_ref_pos2, 1, "omitnan"))./mean(s_ref_pos1, 1, "omitnan");
disp('End of signal position Difference')

end


