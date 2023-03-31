function TC_data = get_tipping_curve_data_miawara(rawSpectra, logFile, calibrationTool)
disp('tipping curve miawara')
% idx hot
idx_hot = find(logFile.Mirror_pos == 0);

% the tipping curve is calculated for the average of a certain frequency
% range
disp('calibrationTool location')

f = load(calibrationTool.channel_freqs);
%disp(f)
N        = calibrationTool.numberOfChannels;
if calibrationTool.read_freq_direct == true
    freq = f(:,2)';
else
    freq     = interp1(f(:,1)',f(:,2)',1:N);  
end
idx_freq = find (freq > 22.135e9 & freq < 22.335e9);

disp('len freq')
disp(length(freq))
disp('len inx freq')
disp(length(idx_hot))

% cycle through hot measurements, the hot measurement is always the first
% measurement

counter = 1;
clear TC_data
%%
for k= 1:length(idx_hot)
    
    disp('starting loop')
    disp(k)

    % find tipping curve measurements within 3min of the hot
    tk =  logFile.time(idx_hot(k));
    
    rainhood = logFile.rainhood_open | logFile.rainhood_open == logFile.RH_closed;

    %disp(logFile.time-tk)
    if k > length(idx_hot)
        tk1 = logFile.time(idx_hot(k+1));
        idx_tmp = find (logFile.Mirror_pos == 2 & logFile.time-tk < 3/60/24 & logFile.time > tk & logFile.time < tk1 & rainhood );
        coldind_tmp = find (logFile.Mirror_pos == 4 & logFile.time-tk < 3/60/24 & logFile.time > tk & logFile.time < tk1 & rainhood );
    else
        idx_tmp = find (logFile.Mirror_pos == 2  & logFile.time-tk < 3/60/24  & rainhood );
        coldind_tmp = find (logFile.Mirror_pos == 4  & logFile.time-tk < 3/60/24 & rainhood );
    end
    

    if (length(idx_tmp)) < 1 || ( length(coldind_tmp) < 1)
        disp('Not enough tipping angles')
        idx_cold    = [];
        idx_tipping = [];
        continue;
    else
        idx_cold    = coldind_tmp(end);
        idx_tipping = idx_tmp(1:end);
    end
    
    % count tipping curve angles
    if length(idx_tipping) < 2
        disp('Not enough tipping angles err2')
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
    
    TC_data{counter}.time = tk;%AB moved to here otherwise can create error

    % select freqency interval line 
    TC_data{counter}.s_tipping = s_tipping (:,idx_freq);
    disp('s_tipping')
    disp(length(TC_data{counter}.s_tipping))
        
    TC_data{counter}.s_hot = s_hot(idx_freq);
    TC_data{counter}.s_cold = s_cold (idx_freq);
    TC_data{counter}.za_tipping = 90 - logFile.Mirror_elevation(idx_tipping);% - el_corr_i;
    TC_data{counter}.za_cold = 90 - logFile.Mirror_elevation(idx_cold);% - el_corr_i;
    TC_data{counter}.frequency = freq(idx_freq);

    % find Thot and Tamb
    disp('finding thot and tamb')
    
    Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7];
    Thot_in     = Thot_in_tmp(:,idx_hot(k));
    Thot_time   = logFile.time(idx_hot(k));
    Tamb_in     = [logFile.meteo.temperature];
    Tamb_time   = [logFile.meteo.time];
    
    [TC_data{counter}.Thot, TC_data{counter}.Tamb] = check_Thot_and_Tamb(Thot_time, Thot_in, Tamb_time, Tamb_in); 
    %clear variables
    
    clear idx_cold idx_tipping
    disp(['counter: ' num2str(counter)])
    %increment counter
    counter = counter + 1;
    
end
disp(['TC len: ' num2str(length(TC_data))])

    

