function [calibratedSpectra] = balancing_calibration_miawara(rawSpectra,logFile,calibrationTool,calType)

N = size(rawSpectra,2);
pol = length(logFile.TC);
m_loop = 0;
calibVersion = calibrationTool.calibrationVersion;
calibrationTool.errorMax = 8;

% read time and mirror position for ref hot and cold and for the whole day

[drift, logFile] = make_drift_miawara(rawSpectra,logFile,calibrationTool);
valInd = ~logFile.isOutlier;

%hot and cold times
t_hot_dis1 = logFile.time(logFile.isHot & logFile.isMirrorDisp1);
t_hot_dis2 = logFile.time(logFile.isHot & logFile.isMirrorDisp2);

t_cold_dis1 = logFile.time(logFile.isColdSky & logFile.isMirrorDisp1 & valInd);
t_cold_dis2 = logFile.time(logFile.isColdSky & logFile.isMirrorDisp2 & valInd);

%hot and cold counts
S_hot1   = rawSpectra(logFile.isHot & logFile.isMirrorDisp1 & valInd,:);
S_hot2   = rawSpectra(logFile.isHot & logFile.isMirrorDisp2 & valInd,:);

S_cold1   = rawSpectra(logFile.isColdSky & logFile.isMirrorDisp1 & valInd,:);
S_cold2   = rawSpectra(logFile.isColdSky & logFile.isMirrorDisp2 & valInd,:);

%mirror cold
mirror_cold_d1 = logFile.Mirror_elevation(logFile.isColdSky & logFile.isMirrorDisp1 & valInd);
mirror_cold_d2 = logFile.Mirror_elevation(logFile.isColdSky & logFile.isMirrorDisp2 & valInd);

%polarisation
j = 1;

if ~isempty(S_hot1) && ~isempty(S_cold1)

    % Check if day before and day after exist
    
    f_pre  = [calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' datestr(datenum(calibrationTool.dateStr)-1,'YYYY_mm_dd') '_' num2str(j) '.txt'];
    f_post = [calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' datestr(datenum(calibrationTool.dateStr)+1,'YYYY_mm_dd') '_' num2str(j) '.txt'];
    
    % load the day before and the day after
    
    if exist(f_pre,'file')
        pre  = readtable(f_pre);
    elseif exist([calibrationTool.rawFileFolder calibrationTool.instrumentName '_' datestr(datenum(calibrationTool.dateStr)-1,'YYYY_mm_dd') '.bin'],'file')    
        write_tau_and_Teff(datenum(calibrationTool.dateStr)-1,calibrationTool.instrumentName);
        pre  = readtable(f_pre);
    else 
        pre.time = logFile.time(1)-1;
        pre.tau  = logFile.TC(j).tau(1);
        pre.Teff = logFile.TC(j).Teff(1);
    end
    
    
    if exist(f_post,'file')
        post  = readtable(f_post);
    elseif exist([calibrationTool.rawFileFolder calibrationTool.instrumentName '_' datestr(datenum(calibrationTool.dateStr)+1,'YYYY_mm_dd') '.bin'],'file')    
        write_tau_and_Teff(datenum(calibrationTool.dateStr)+1,calibrationTool.instrumentName);
        post  = readtable(f_post);
    else 
        post.time = logFile.time(1)+1;
        post.tau  = logFile.TC(j).tau(end);
        post.Teff = logFile.TC(j).Teff(end);
        %check Tdiff ?? min?
    end
    
    
    %check if tau is 0
    
    t    = [pre.time' logFile.TC(j).time post.time'];
    tau  = [pre.tau'  logFile.TC(j).tau post.tau'];
    Teff = [pre.Teff' logFile.TC(j).Teff post.Teff'];

    %t_corr   = t(tau ~= 0 ); %lots of nans sometimes included 
    %tau_corr = tau(tau ~= 0);%tau without extra condition
    
    t_corr   = t((tau ~= 0) & ~isnan(tau)); %lots of nans sometimes included 
    tau_corr = tau((tau ~= 0) & ~isnan(tau));%tau without extra condition
    
       
    %Time periods over which to make the calibration
    t0 = datenum(calibrationTool.dateStr,'yyyy_mm_dd');
    dt_int = calibrationTool.calibrationTime / 60 / 24 ;
    
    if ~isempty(rawSpectra)
      % loop over line measurements
    
      for m = 1:1/dt_int
        %change loop time period variables  
        m_loop = m_loop+1;
        t1 = t0+(m-1)*dt_int;
        t2 = t0+m*dt_int;
    
        disp(['run calibration for ' datestr(t1,'yyyy-mm-dd HH:MM') ' - ' datestr(t2,'yyyy-mm-dd HH:MM')])
    
        %time condition for selecting data
        isTint = logFile.time >= t1 & logFile.time < t2;
    
        %time of valid line and ref measurement 
        t_line_dis1 = logFile.time(isTint & logFile.isLine & logFile.isMirrorDisp1 & valInd);
        t_line_dis2 = logFile.time(isTint & logFile.isLine & logFile.isMirrorDisp2 & valInd);
    
        t_ref_dis1 = logFile.time(isTint & logFile.isRef & logFile.isMirrorDisp1 & valInd);
        t_ref_dis2 = logFile.time(isTint & logFile.isRef & logFile.isMirrorDisp2 & valInd);
    
        %mirror line and ref 
        mirror_ref_d1  = logFile.Mirror_elevation(isTint &logFile.isRef & logFile.isMirrorDisp1 & valInd);
        mirror_ref_d2  = logFile.Mirror_elevation(isTint &logFile.isRef & logFile.isMirrorDisp2 & valInd);
    
        mirror_line_d1 = logFile.Mirror_elevation(isTint & logFile.isLine & logFile.isMirrorDisp1 & valInd);
        mirror_line_d2 = logFile.Mirror_elevation(isTint & logFile.isLine & logFile.isMirrorDisp2 & valInd);
    
        %Counts line and Ref
        S_line1   = rawSpectra(isTint & logFile.isLine & logFile.isMirrorDisp1 & valInd,:);
        S_line2   = rawSpectra(isTint & logFile.isLine & logFile.isMirrorDisp2 & valInd,:);
    
        S_ref1   = rawSpectra(isTint & logFile.isRef & logFile.isMirrorDisp1 & valInd,:);
        S_ref2   = rawSpectra(isTint & logFile.isRef & logFile.isMirrorDisp2 & valInd,:);
        
        %interpolation of optical depth
        TAU = interp1(t_corr,tau_corr,t_line_dis1);
        TEFF = interp1(t,Teff,t_line_dis1);

        for k = 1:size(S_line1,1)
    
            clear data
    
            % find closest to displacement1 sky measurement for displacement2 
            % measurement, the index closest to the displacement1 measurement 
            % is found for view i.e. hotdis2 closest to hotdis1, refdis2 closest
            %to refdis1
            [dt_line2,  idx_line2(k)]  = min(abs( t_line_dis1(k) - t_line_dis2));
    
            [dt_ref1,  idx_ref1(k)]  = min(abs( t_line_dis1(k) - t_ref_dis1));
            [dt_ref2,  idx_ref2(k)]  = min(abs( t_ref_dis1(idx_ref1(k)) - t_ref_dis2));
    
            [dt_hot1,  idx_hot1(k)]  = min(abs( t_line_dis1(k) - t_hot_dis1));
            [dt_hot2,  idx_hot2(k)]  = min(abs( t_hot_dis1(idx_hot1(k)) - t_hot_dis2));
    
            [dt_cold1, idx_cold1(k)] = min(abs(t_line_dis1(k) - t_cold_dis1));
            [dt_cold2, idx_cold2(k)] = min(abs(t_cold_dis1(idx_cold1(k)) - t_cold_dis2));
    
            %conditions to ensure measurements are matching
            nonempty_cond = ~isempty(dt_hot1) && ~isempty(dt_ref1) && ~isempty(dt_cold1) ...
            && ~isempty(dt_hot1) && ~isempty(dt_ref1) && ~isempty(dt_cold1) ;
            time_cond = (dt_ref1*24*60<2*60) && (dt_cold1*24*60<40*60) && (dt_hot1*24*60<40*60) ;
            mirror_cond = (mirror_line_d1(k) == mirror_line_d2(idx_line2(k))) && ...
                (mirror_cold_d1(idx_cold1(k)) == mirror_cold_d2(idx_cold2(k)));
    
            % check delta t
            if nonempty_cond && time_cond && mirror_cond
    
                % find Thot and Tamb
                if strcmp( calibrationTool.instrumentName,'MIAWARA' )
                    Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7];
                else
                    Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7;logFile.THot8;logFile.THot9];
                end
    
    
                Thot_in     =  Thot_in_tmp(:,idx_hot1(k))';
                Thot_time   = logFile.time(idx_hot1(k));
                Tamb_in     = [logFile.meteo.temperature];
                Tamb_time   = [logFile.meteo.time];
    
                [Thot(k), ~] = check_Thot_and_Tamb(Thot_time, Thot_in, Tamb_time, Tamb_in); 
    
    
                % write data into struct
    
                data.S_line = mean( [S_line1(k,:); S_line2(idx_line2(k),:)],1 ) ;
                data.S_ref  = mean( [S_ref1(idx_ref1(k),:); S_ref2(idx_ref2(k),:)],1 ) ;
                data.S_hot  = mean( [S_hot1(idx_hot1(k),:); S_hot2(idx_hot2(k),:)],1 ) ;
                data.S_cold = mean( [S_cold1(idx_cold1(k),:); S_cold2(idx_cold2(k),:)],1 ) ;
                data.Thot   = Thot(k);
                data.Tref   = Thot(k);
                data.Teff   = TEFF(k);
                data.T0     = 2.7;
                data.tau    = TAU(k);
                data.mirror_el_line = mirror_line_d1(k);
                data.mirror_el_ref  = mirror_ref_d1(idx_ref1(k));
                data.mirror_el_cold = mirror_cold_d1(idx_cold1(k));
                if isfield(calibrationTool,'errorMax')
                  data.errorMax = calibrationTool.errorMax;
                end
    
                % counts2kelvin
                if calibrationTool.DoColdCal == 1
                    [out.Tb(k,:),out.T_COLD(k,:),out.T_rec(k,:),out.A(k,:),out.a(k,:),out.delta_Tb(k,:)] = counts2kelvin (data);    % offset
                    Tb_cold_calibrated(k,:) = coldcounts2kelvin(data);
                else
                    [out.Tb(k,:),out.T_COLD(k,:),out.T_rec(k,:),out.A(k,:),out.a(k,:),out.delta_Tb(k,:)] = counts2kelvin (data);    % offset
                end
               
            else
                disp('Measurement mismatch: skipping this one')
              
            end
    
        end
        if calibrationTool.DoColdCal
            Tb_cold_calibrated_mean = mean(Tb_cold_calibrated,1,'omitnan');
            save_tb_coldcal_miawara(calibrationTool,  Tb_cold_calibrated_mean, num2str(m))
        end
        calibratedSpectra(m_loop).Tbnancount = nnz(isnan(out.Tb(:,1)));
        calibratedSpectra(m_loop).Tb = mean(out.Tb, 'omitnan');
        calibratedSpectra(m_loop).stdTb = std(out.Tb, 'omitnan'); 
        calibratedSpectra(m_loop).meanStdTb = mean(std(out.Tb, 'omitnan'), 'omitnan');
        calibratedSpectra(m_loop).TSys = mean(mean(out.T_rec, 'omitnan'), 'omitnan');
        calibratedSpectra(m_loop).stdTSys = mean(std(out.T_rec, 'omitnan'), 'omitnan');
        calibratedSpectra(m_loop).tau = mean(TAU, 'omitnan');
        calibratedSpectra(m_loop).A   = mean(out.A, 'omitnan');     % A:      the equivalent transmission of the reference absorber
        calibratedSpectra(m_loop).a   = mean(out.a, 'omitnan');     % a:      the correction coefficient for the troposphere and ref absorber
        calibratedSpectra(m_loop).TCold = mean(out.T_COLD, 'omitnan');
                             
        %calibratedSpectra(m_loop).antennaIndCleanAngle = isLineTimeDir;
        calibratedSpectra(m_loop).refInd  = idx_ref1;
        calibratedSpectra(m_loop).hotInd  = idx_hot1;
        calibratedSpectra(m_loop).coldInd = idx_cold1;
    
        clear idx_ref idx_hot idx_cold
    
        calibratedSpectra(m_loop).THot  = mean(Thot, 'omitnan');
        calibratedSpectra(m_loop).TCold = mean(out.T_COLD, 'omitnan');
        calibratedSpectra(m_loop).theoreticalStartTime = (m-1)*dt_int *24; % in hours
        calibratedSpectra(m_loop).calibrationTime      = calibrationTool.calibrationTime;
        calibratedSpectra(m_loop).calibrationVersion   = calibVersion;
    
        clear out
      end      
  else
    disp('No data for this time interval!')
  end
end

end
