function [calibratedSpectra] = balancing_calibration_generic(rawSpectra,logFile,calibrationTool,calType)
%==========================================================================
% NAME          | balancing_calibration_generic
% TYPE          | function
% AUTHOR(S)     | Franzisca Schranz
% CREATION      | 06.2020
%               |
% ABSTRACT      | ...
%               |
%               |
% ARGUMENTS     | INPUTS:  1. rawSpectra
%               |          2. logFile
%               |          3. calibrationTool
%               |            - calibrationVersion
%               |            - ...
%               |          4. calType
%               |
%               | OUTPUTS: 1. calibratedSpectra
%               |
% COMMENTS      |
%               |
%==========================================================================

% Calibration version
calibVersion = calibrationTool.calibrationVersion ;

% observation direction

if strcmp(calType,'left')
    isDir = logFile.dir == 1;
elseif strcmp(calType,'right')
    isDir = logFile.dir == 2;
else
    isDir = ones(size(logFile.Year));
end

%%
% 
N = size(rawSpectra,2);
pol = length(logFile.TC);
  
% read time and mirror position for ref hot and cold and for the whole day

t_ref   = logFile.time(logFile.isRef & isDir);  % select ref with specific direction left ==1 or right ==2
t_hot   = logFile.time(logFile.isHot);
t_cold  = logFile.time(logFile.isColdSky);

mirror_ref  = logFile.Mirror_elevation(logFile.isRef & isDir);
mirror_cold = logFile.Mirror_elevation(logFile.isColdSky);


for j = 1:length(logFile.TC) % loop over polarisations
    
    
    % load remaining data from day before
    
    %%%%[logFile,rawSpectra] = calibrationTool.read_level0(calibrationTool);
    
    
    
    % define indices for polarisation (if polarized measurement: first half of spectrum:pol1; second half of spectrum:pol2)
    
    a = (j-1)*N/pol + 1;
    b = j*N/pol;
    
    
    % read ref, hot and cold spectra for 1 polarisation and the whole day
        
    S_ref   = rawSpectra(logFile.isRef & isDir,a:b);
    S_hot   = rawSpectra(logFile.isHot,a:b);
    S_cold  = rawSpectra(logFile.isColdSky,a:b);

        
    % % Prepare opacity and Teff
    
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
        check Tdiff ?? min?
    end
    
    
    %check if tau is 0
    
    t    = [pre.time' logFile.TC(j).time post.time'];
    tau  = [pre.tau'  logFile.TC(j).tau post.tau'];
    Teff = [pre.Teff' logFile.TC(j).Teff post.Teff'];
    
    t_corr   = t(tau ~= 0);
    tau_corr = tau(tau ~= 0);
       
           
    % prepare loop over integration time intervals 
    
    t0 = datenum(calibrationTool.dateStr,'yyyy_mm_dd');
    dt_int = calibrationTool.calibrationTime / 60 / 24 ;
    
    
    % loop over integration time intervals
    
    for m = 1:1/dt_int
        
        t1 = t0+(m-1)*dt_int;
        t2 = t0+m*dt_int;
        
        % Read line spectra within time interval

        isInTimeInterval = logFile.time >= t1 & logFile.time < t2;
        isLineTimeDir    = logFile.isLine & isInTimeInterval & isDir;
        
        t_line      = logFile.time(isLineTimeDir);
        mirror_line = logFile.Mirror_elevation(isLineTimeDir);
        S_line      = rawSpectra(isLineTimeDir,a:b);
        
        TAU = interp1(t_corr,tau_corr,t_line);
        TEFF = interp1(t,Teff,t_line);

        
        % loop over line measurements
        
        for k = 1:size(S_line,1)
    
            clear data

            % find closest

            [dt_ref,  idx_ref(k)]  = min(abs( t_line(k) - t_ref));
            [dt_hot,  idx_hot(k)]  = min(abs( t_line(k) - t_hot));
            [dt_cold, idx_cold(k)] = min(abs( t_line(k) - t_cold));


            % check delta t

            if ~isempty(dt_hot) && ~isempty(dt_ref) && ~isempty(dt_cold)
                if (dt_ref*24*60<2*60) && (dt_cold*24*60<40*60) && (dt_hot*24*60<40*60)



                    % find Thot and Tamb

                    Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7;logFile.THot8;logFile.THot9];
                    Thot_in     =  Thot_in_tmp(:,idx_hot(k))';
                    Thot_time   = logFile.time(idx_hot(k));
                    Tamb_in     = logFile.meteo.temperature;
                    Tamb_time   = logFile.meteo.time;

                    [Thot(k), ~] = check_Thot_and_Tamb(Thot_time, Thot_in, Tamb_time, Tamb_in); 
                    

                    % write data into struct

                    data.S_line = S_line(k,:);
                    data.S_ref  = S_ref(idx_ref(k),:);
                    data.S_hot  = S_hot(idx_hot(k),:);
                    data.S_cold = S_cold(idx_cold(k),:);
                    data.Thot   = Thot(k);
                    data.Tref   = Thot(k);
                    data.Teff   = TEFF(k);
                    data.T0     = 2.7;
                    data.tau    = TAU(k);
                    data.mirror_el_line = mirror_line(k);
                    data.mirror_el_ref  = mirror_ref(idx_ref(k));
                    data.mirror_el_cold = mirror_cold(idx_cold(k));


                    % counts2kelvin
                    %[Tb,T_COLD,T_rec,error,A,a,delta_Tb]=counts2kelvin_v6(S_LINE,S_REF,S_HOT,S_COLD,T_HOT,T_REF,T_EFF,T0,TAU,mirror_elevation_line,mirror_elevation_ref,mirror_elevation_cold,error);
                    [out.Tb(k,:),out.T_COLD(k,:),out.T_rec(k,:),out.A(k,:),out.a(k,:),out.delta_Tb(k,:)] = counts2kelvin (data);    % offset



                end
            end
            
        end
        
        % write measurements of day after to file

        isInTimeInterval = logFile.time >= t0+1;
        
        f = [calibrationTool.file '_missing'];
        
        M = rawSpectra(isInTimeInterval,:)';
        fid = fopen([f '.bin'], 'w');
        fwrite(fid,M(:));
        fclose(fid);
        
        writecell(logFile.header',[f '.txt'],'Delimiter',';')
        dlmwrite([f '.txt'],logFile.x(:,isInTimeInterval)','delimiter',';','-append');
        
        % save averaged spectra
        
        m_loop = m + (j-1)/dt_int;
        
        calibratedSpectra(m_loop).antennaIndCleanAngle = isLineTimeDir;
        calibratedSpectra(m_loop).refInd  = idx_ref;
        calibratedSpectra(m_loop).hotInd  = idx_hot;
        calibratedSpectra(m_loop).coldInd = idx_cold;
        
        clear idx_ref idx_hot idx_cold
        
        calibratedSpectra(m_loop).THot  = nanmean(Thot);
        calibratedSpectra(m_loop).TCold = nanmean(out.T_COLD);
        
        
        calibratedSpectra(m_loop).Tb = nanmean(out.Tb);
        calibratedSpectra(m_loop).TSys = nanmean(nanmean(out.T_rec));
        calibratedSpectra(m_loop).stdTSys = nanmean(std(out.T_rec));
        
%         calibratedSpectra(m).Yspectral=nanmean(S_hot(idx_hot))./nanmean(S_cold(idx_cold));
%         calibratedSpectra(m).TN=(calibratedSpectra(m).THot - calibratedSpectra(m).Yspectral*calibrationTool.TCold)./ (calibratedSpectra(m).Yspectral -1);
        calibratedSpectra(m_loop).pol = j;
        calibratedSpectra(m_loop).dir = logFile.dir(find(isDir ==1,1,'first'));
        calibratedSpectra(m_loop).theoreticalStartTime = (m-1)*dt_int *24; % in hours
        calibratedSpectra(m_loop).calibrationTime=calibrationTool.calibrationTime;
        calibratedSpectra(m_loop).calibrationVersion = calibVersion;
        clear out
    end
    
    
    % mean ?min
%left right

end



    