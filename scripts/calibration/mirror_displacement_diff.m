function mirror_displacement_diff(rawSpectra,logFile,calibrationTool,calType)
%==========================================================================
% NAME          | mirror_displacement_diff
% TYPE          | function
% AUTHOR(S)     | Alistair Bell adapted from balancing calibration by Franzisca Schranz
% CREATION      | 11.2022
%               |
% ARGUMENTS     | INPUTS:  1. ToEdit
%               |          
%               |
%               | OUTPUTS: 1. ToEdit
%               |
% COMMENTS      |To Add
%               |
%==========================================================================
%physical constants
disp('starting mirror displacement diff')
RE = 6371000;


% Calibration version
%calibVersion = calibrationTool.calibrationVersion;

% base index for data storage
m_loop = 0;

% observation direction


isDir = ones(size(logFile.Year));

%%
N = size(rawSpectra,2);
pol = length(logFile.TC);

% read time and mirror position for ref hot and cold and for the whole day

t_ref   = logFile.time(logFile.isRef & isDir);  % select ref with specific direction left ==1 or right ==2
t_hot   = logFile.time(logFile.isHot);
t_cold  = logFile.time(logFile.isColdSky);

mirror_ref  = logFile.Mirror_elevation(logFile.isRef & isDir);

mirror_cold = logFile.Mirror_elevation(logFile.isColdSky);
       
% read ref, hot and cold spectra for 1 polarisation and the whole day

S_ref_disp1   = rawSpectra(logFile.isRef & logFile.isMirrorDisp1, :);
S_hot_disp1   = rawSpectra(logFile.isHot & logFile.isMirrorDisp1, :);
S_cold_disp1  = rawSpectra(logFile.isColdSky & logFile.isMirrorDisp1, :);

S_ref_disp2   = rawSpectra(logFile.isRef & logFile.isMirrorDisp2, :);
S_hot_disp2   = rawSpectra(logFile.isHot & logFile.isMirrorDisp2,:);
S_cold_disp2  = rawSpectra(logFile.isColdSky & logFile.isMirrorDisp2, :);

t_ref_disp1   = logFile.time(logFile.isRef & isDir & logFile.isMirrorDisp1);  
t_hot_disp1   = logFile.time(logFile.isHot & logFile.isMirrorDisp1);
t_cold_disp1  = logFile.time(logFile.isColdSky & logFile.isMirrorDisp1);

t_ref_disp2   = logFile.time(logFile.isRef & isDir & logFile.isMirrorDisp2);  
t_hot_disp2   = logFile.time(logFile.isHot & logFile.isMirrorDisp2);
t_cold_disp2  = logFile.time(logFile.isColdSky & logFile.isMirrorDisp2);

mirror_ref_disp1  = logFile.Mirror_elevation(logFile.isRef &  logFile.isMirrorDisp1);
mirror_ref_disp2  = logFile.Mirror_elevation(logFile.isRef &  logFile.isMirrorDisp2);

mirror_cold = logFile.Mirror_elevation(logFile.isColdSky &  logFile.isMirrorDisp1);
mirror_cold = logFile.Mirror_elevation(logFile.isColdSky &  logFile.isMirrorDisp2);

%%
if ~isempty(S_ref_disp1) & ~isempty(S_hot_disp1) & ~isempty(S_hot_disp1) & ...
~isempty(S_ref_disp2) & ~isempty(S_hot_disp2) & ~isempty(S_hot_disp2)
    
    % Prepare opacity and Teff

    % Check if day before and day after exist

    % load the day before and the day after

        t    = logFile.TC.time;
        tau  = logFile.TC.tau;
        Teff = logFile.TC.Teff;
        
        %check if tau is 0

        t_corr   = t(tau ~= 0);
        tau_corr = tau(tau ~= 0);

        % prepare loop over integration time intervals 

        t0 = datenum(calibrationTool.dateStr,'yyyy_mm_dd');
        dt_int = calibrationTool.calibrationTime / 60 / 24 ;


        % loop over integration time intervals

        for m = 1:1/dt_int

            t1 = t0+(m-1)*dt_int; %define time bounds
            t2 = t0+m*dt_int;

            disp(['run calibration for ' datestr(t1,'yyyy-mm-dd HH:MM') ' - ' datestr(t2,'yyyy-mm-dd HH:MM')])

            % Read line spectra within time interval

            isInTimeInterval = logFile.time >= t1 & logFile.time < t2;
            isLineTimeDir    = logFile.isLine & isInTimeInterval & isDir;

            t_line      = logFile.time(isLineTimeDir);
            mirror_line = logFile.Mirror_elevation(isLineTimeDir);
            S_line      = rawSpectra(isLineTimeDir,:);

            %select reference measurements that are in the initial time
            %interval
            t_ref_disp1_intime = logFile.time(logFile.isRef & isInTimeInterval & logFile.isMirrorDisp1); 
            s_ref_disp1_intime = rawSpectra( logFile.isRef & isInTimeInterval & logFile.isMirrorDisp1,:);

            %mirror_line = logFile.Mirror_elevation(isLineTimeDir);
          
            mirror_ref_el_disp1 = logFile.Mirror_elevation(logFile.isRef & isInTimeInterval & logFile.isMirrorDisp1);

            %stop from interpolating for tau before tipping curve has been
            %made
                    % Check if day before and day after exist
            j = 1;
            f_pre  = [calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' datestr(datenum(calibrationTool.dateStr)-1,'YYYY_mm_dd') '_' num2str(1) '.txt'];
            f_post = [calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' datestr(datenum(calibrationTool.dateStr)+1,'YYYY_mm_dd') '_' num2str(1) '.txt'];
    
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
    
            t_corr   = t(tau ~= 0);
            tau_corr = tau(tau ~= 0);

            TAU = interp1(t_corr,tau_corr,t_ref_disp1_intime);
            TEFF = interp1(t,Teff,t_ref_disp1_intime);
            
            if ~isempty(s_ref_disp1_intime)
                % loop over line measurements
                m_loop = m_loop+1;
                
                for k = 1:size(s_ref_disp1_intime,1) %for each line measurement, nearest hot, cold and reference are found for each spectral reading

                    clear data

                    % find closest                   
                    [dt_ref_disp2,  idx_ref_disp2(k)]  = min(abs( t_ref_disp1_intime(k) - t_ref_disp2));
                    
                    [dt_hot_disp1,  idx_hot_disp1(k)]  = min(abs( t_ref_disp1_intime(k) - t_hot_disp1));
                    [dt_hot_disp2,  idx_hot_disp2(k)]  = min(abs( t_ref_disp1_intime(k) - t_hot_disp2));
    
                    [dt_cold_disp1,  idx_cold_disp1(k)]  = min(abs( t_ref_disp1_intime(k) - t_hot_disp1));
                    [dt_cold_disp2,  idx_cold_disp2(k)]  = min(abs( t_ref_disp1_intime(k) - t_cold_disp1));

                    disp('hot dis pos. 1 idx; Ref dis pos. 2 idx ')
                    disp([idx_hot_disp1(k),  idx_ref_disp2(k)])
                
                    % check delta t

                    if ~isempty(dt_ref_disp2) && ~isempty(dt_hot_disp1) && ~isempty(dt_cold_disp1)
                        if (dt_ref_disp2*24*60<2*60) && (dt_hot_disp1*24*60<40*60) && (dt_cold_disp1*24*60<40*60)

                            % find Thot and Tamb
                            if strcmp( calibrationTool.instrumentName,'MIAWARA' )
                                Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7];
                            else
                                Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7;logFile.THot8;logFile.THot9];
                            end


                            Thot_disp1_in = Thot_in_tmp(:,idx_hot_disp1(k))';
                            Thot_disp1_time   = logFile.time(idx_hot_disp1(k));
                                
                            Thot_disp2_in = Thot_in_tmp(:,idx_hot_disp2(k))';
                            Thot_disp2_time   = logFile.time(idx_hot_disp2(k));

                            Tamb_in     = [logFile.meteo.temperature];
                            Tamb_time   = [logFile.meteo.time];

                            [Thot_disp1(k), ~] = check_Thot_and_Tamb(Thot_disp1_time, Thot_disp1_in, Tamb_time, Tamb_in); 
                            [Thot_disp2(k), ~] = check_Thot_and_Tamb(Thot_disp2_time, Thot_disp2_in, Tamb_time, Tamb_in); 

                            % write data into struct

                            data.S_line = S_line(k,:);

                            data.S_ref_disp1  = s_ref_disp1_intime(k,:);
                            data.S_ref_disp2  = S_ref_disp2(idx_ref_disp2(k),:);
                            
                            data.S_hot_disp1  = S_hot_disp1(idx_hot_disp1(k),:);
                            data.S_hot_disp2  = S_hot_disp2(idx_hot_disp2(k),:);

                            data.S_cold_disp1 = S_cold_disp1(idx_cold_disp1(k),:);
                            data.S_cold_disp2= S_cold_disp2(idx_cold_disp2(k),:);

                            data.Thot_disp1   = Thot_disp1(k);
                            data.Thot_disp2   = Thot_disp2(k);

                            data.Tref_disp1   = Thot_disp1(idx_hot_disp1(k));
                            data.Tref_disp2   = Thot_disp1(idx_hot_disp1(k));
                            
                            data.Teff   = TEFF(k);
                            data.T0     = 2.7;
                            data.tau    = TAU(k);

                            data.mirror_el_line = mirror_line(k);
                            data.mirror_el_ref  = mirror_ref(k);

                            data.mirror_el_cold_disp1 = mirror_cold(idx_cold_disp1(k));
                            data.mirror_el_cold_disp2 = mirror_cold(idx_cold_disp2(k));
                            
                            [Tb_ref_diff(k,:),Tb_hot_diff(k,:), Tb_cold_diff(k,:)] = displacementcounts2kelvin(data);
                            end
                    end
                end
                
            end
            Tb_ref_diff_mean = mean(Tb_ref_diff, 1, 'omitnan');
            Tb_hot_diff_mean = mean(Tb_hot_diff, 1, 'omitnan');
            Tb_cold_diff_mean = mean(Tb_cold_diff, 1, 'omitnan');

            save_difference_data_miawara(calibrationTool, Tb_hot_diff_mean, Tb_ref_diff_mean, Tb_cold_diff_mean, int2str(m));
        end
    else
        disp('No data!')
    end
end
    
    


    