
function TC_data = get_tipping_curve_data_gromos(rawSpectra, logFile, calibrationTool)


% idx hot

idx_hot = find(logFile.Mirror_pos == 1);


% the tipping curve is calculated for the average of a certain frequency
% range

f = load(calibrationTool.channel_freqs);
N        = calibrationTool.numberOfChannels;
freq     = interp1(f(:,1)',f(:,2)',1:N/2);
idx_freq = find (freq > 22.135e9 & freq < 22.335e9);


% cycle through hot measurements, the hot measurement is always the first
% measurement

counter = 1;
clear TC_data
%%
for k= 1:length(idx_hot)
    
    % find tipping curve measurements within 3min of the hot
    tk =  logFile.time(idx_hot(k));
    
    TC_data{counter}.time = tk;
    
    rainhood = logFile.rainhood_open | logFile.rainhood_open == logFile.RH_closed;
    if k > length(idx_hot)
        tk1 = logFile.time(idx_hot(k+1));
        idx_tmp = find (logFile.Mirror_pos == 6 & logFile.time-tk < 3/60/24 & logFile.time > tk & logFile.time < tk1 & rainhood ); 
    else
        idx_tmp = find (logFile.Mirror_pos == 6 & logFile.time-tk < 3/60/24 & logFile.time > tk & rainhood ); 
    end
    
    idx_cold    = idx_tmp(end);
    idx_tipping = idx_tmp(1:end-1);

    
    % count tipping curve angles
    
    if length(idx_tipping) < 2
        disp('Not enough tipping angles')
        continue;
    end
    
    
    % get lines
        
    s_tipping = rawSpectra (idx_tipping,:);
    s_hot     = rawSpectra (idx_hot(k),:);
    s_cold    = rawSpectra (idx_cold,:);
    
    
    % check lines
    
    ind = find (s_tipping == 0);
    s_tipping (ind) = NaN;
       
    ind = find (s_hot == 0);
    s_hot (ind) = NaN;
       
    ind = find (s_cold == 0);
    s_cold (ind) = NaN;
    
    
    % select freqency interval line 
    
    TC_data{counter}.s_tipping_pol1 = s_tipping (:,idx_freq);
    TC_data{counter}.s_tipping_pol2 = s_tipping (:,idx_freq + N/2);
    
    TC_data{counter}.s_hot_pol1 = s_hot (idx_freq);
    TC_data{counter}.s_hot_pol2 = s_hot (idx_freq + N/2);
    
    TC_data{counter}.s_cold_pol1 = s_cold (idx_freq);
    TC_data{counter}.s_cold_pol2 = s_cold (idx_freq + N/2);
    
    
    % elevation correction time dependent
    
    load(calibrationTool.elcorr_file)
    el_corr_i = interp1 (elcorr(:,1), elcorr(:,2), logFile.time(k));
    
    TC_data{counter}.za_tipping = 90 - logFile.Mirror_elevation(idx_tipping) - el_corr_i;
    TC_data{counter}.za_cold    = 90 - logFile.Mirror_elevation(idx_cold) - el_corr_i;
    
    
    % find Thot and Tamb
    
    Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7;logFile.THot8;logFile.THot9];
    Thot_in     =  Thot_in_tmp(:,idx_hot(k))';
    Thot_time   = logFile.time(idx_hot(k));
    Tamb_in     = logFile.meteo.temperature;
    Tamb_time   = logFile.meteo.time;

    [TC_data{counter}.Thot, TC_data{counter}.Tamb] = check_Thot_and_Tamb(Thot_time, Thot_in, Tamb_time, Tamb_in); 
        
    
    % clear variables
    
    clear idx_cold idx_tipping
    
    % increment counter
    counter = counter + 1;
    
end
    

