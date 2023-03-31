function [o, calibratedSpectra] = run_LN2_calibration_miawara(rawSpectra,logFile,calibrationTool,~)

N = size(rawSpectra,2);
m_loop = 0;
calibVersion = calibrationTool.calibrationVersion;
calibrationTool.errorMax = 8;
T_LN2 = 80;
o = [];
nan_rep = ones(1,N)*nan;


%hot, cold and LN2 times
t_hot_dis1 = logFile.time(logFile.isHot & logFile.isMirrorDisp1);
t_hot_dis2 = logFile.time(logFile.isHot & logFile.isMirrorDisp2);

t_LN2_dis1 = logFile.time(logFile.isColdLN2 & logFile.isMirrorDisp1);
t_LN2_dis2 = logFile.time(logFile.isColdLN2 & logFile.isMirrorDisp2);

%hot and cold counts
S_hot1   = rawSpectra(logFile.isHot & logFile.isMirrorDisp1,:);
S_hot2   = rawSpectra(logFile.isHot & logFile.isMirrorDisp2,:);

S_LN2_1   = rawSpectra(logFile.isColdLN2 & logFile.isMirrorDisp1,:);
S_LN2_2   = rawSpectra(logFile.isColdLN2 & logFile.isMirrorDisp2,:);

%mirror cold
mirror_cold_d1 = logFile.Mirror_elevation(logFile.isColdSky & logFile.isMirrorDisp1);
mirror_cold_d2 = logFile.Mirror_elevation(logFile.isColdSky & logFile.isMirrorDisp2);

%polarisation
j = 1;

if ~isempty(S_hot1) && ~isempty(S_LN2_1)
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
        t_line_dis1 = logFile.time(isTint & logFile.isLine & logFile.isMirrorDisp1 );
        t_line_dis2 = logFile.time(isTint & logFile.isLine & logFile.isMirrorDisp2);
    
        t_ref_dis1 = logFile.time(isTint & logFile.isRef & logFile.isMirrorDisp1 );
        t_ref_dis2 = logFile.time(isTint & logFile.isRef & logFile.isMirrorDisp2);

        t_cold_dis1 = logFile.time(isTint & logFile.isColdSky & logFile.isMirrorDisp1);
        t_cold_dis2 = logFile.time(isTint & logFile.isColdSky & logFile.isMirrorDisp2);

        %mirror line and ref 
        mirror_ref_d1  = logFile.Mirror_elevation(isTint &logFile.isRef & logFile.isMirrorDisp1 );
        mirror_ref_d2  = logFile.Mirror_elevation(isTint &logFile.isRef & logFile.isMirrorDisp2);
    
        mirror_line_d1 = logFile.Mirror_elevation(isTint & logFile.isLine & logFile.isMirrorDisp1);
        mirror_line_d2 = logFile.Mirror_elevation(isTint & logFile.isLine & logFile.isMirrorDisp2);
    
        %Counts line and Ref
        S_line1   = rawSpectra(isTint & logFile.isLine & logFile.isMirrorDisp1,:);
        S_line2   = rawSpectra(isTint & logFile.isLine & logFile.isMirrorDisp2,:);
    
        S_ref1   = rawSpectra(isTint & logFile.isRef & logFile.isMirrorDisp1,:);
        S_ref2   = rawSpectra(isTint & logFile.isRef & logFile.isMirrorDisp2 ,:);

        S_cold1   = rawSpectra(isTint & logFile.isColdSky & logFile.isMirrorDisp1,:);
        S_cold2   = rawSpectra(isTint & logFile.isColdSky & logFile.isMirrorDisp2,:);
        

        nonempty_cond_hcLN2 =  ~isempty(t_hot_dis1) && ~isempty(t_hot_dis2) ...
                && ~isempty(t_cold_dis1) && ~isempty(t_cold_dis2) && ...
        ~isempty(t_LN2_dis1) && ~isempty(t_LN2_dis2);

        nonempty_line =  ~isempty(t_line_dis1) && ~isempty(t_line_dis2);
        nonempty_ref =  ~isempty(t_ref_dis1) && ~isempty(t_ref_dis2);

        clear data


        T_LINE = zeros(size(S_cold1,1),N);
        T_REF = zeros(size(S_cold1,1),N);
        T_COLD = zeros(size(S_cold1,1),N);
        T_rec = zeros(size(S_cold1,1),N);

        if nonempty_cond_hcLN2
            for k = 1:size(S_cold1,1)
                disp(k)        
                % find closest to displacement1 sky measurement for displacement2 
                % measurement, the index closest to the displacement1 measurement 
                % is found for view i.e. hotdis2 closest to hotdis1, refdis2 closest
                %to refdis1
                if nonempty_line
                  [dt_line1,  idx_line1(k)]  = min(abs( t_cold_dis1(k) - t_line_dis1));
                  [dt_line2,  idx_line2(k)]  = min(abs(t_line_dis1(idx_line1(k)) - t_line_dis2));
                end

                if nonempty_ref
                  [dt_ref1,  idx_ref1(k)]  = min(abs( t_cold_dis1(k) - t_ref_dis1));
                  [dt_ref2,  idx_ref2(k)]  = min(abs( t_ref_dis1(idx_ref1(k)) - t_ref_dis2));
                end

                [dt_hot1,  idx_hot1(k)]  = min(abs( t_cold_dis1(k) - t_hot_dis1));
                [dt_hot2,  idx_hot2(k)]  = min(abs( t_hot_dis1(idx_hot1(k)) - t_hot_dis2));
        
                [dt_cold1, idx_cold1(k)] = min(abs(t_cold_dis1(k) - t_cold_dis1));
                [dt_cold2, idx_cold2(k)] = min(abs(t_cold_dis1(idx_cold1(k)) - t_cold_dis2));

                [dt_LN2_1,  idx_LN2_1(k)]  = min(abs( t_cold_dis1(k) - t_LN2_dis1));
                [dt_LN2_2,  idx_LN2_2(k)]  = min(abs( t_LN2_dis1(idx_LN2_1(k)) - t_LN2_dis2));
        
                %conditions to ensure measurements are matching
                nonempty_cond = ~isempty(dt_hot1) && ~isempty(dt_hot2) && ~isempty(dt_LN2_1) ...
                && ~isempty(dt_LN2_2);

                nonempty_line =  exist('dt_line1', 'var')  && ~isempty(dt_line1);
                nonempty_ref = exist('dt_ref1', 'var')  && ~isempty(dt_ref1) ;
                nonempty_hot =  exist('dt_hot1', 'var')  &&~isempty(dt_hot1);
                nonempty_cold = exist('dt_cold1', 'var')  && ~isempty(dt_cold1);

                     
                % check delta t
                if nonempty_cond
        
                    % find Thot and Tamb
                    if strcmp( calibrationTool.instrumentName,'MIAWARA' )
                        Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7];
                    else
                        Thot_in_tmp = [logFile.THot0;logFile.THot1;logFile.THot2;logFile.THot3;logFile.THot4;logFile.THot5;logFile.THot6;logFile.THot7;logFile.THot8;logFile.THot9];
                    end
                
                    %Find temperature of THOT
                    Thot_in     =  Thot_in_tmp(:,idx_hot1(k))';
                    Thot_time   = logFile.time(idx_hot1(k));
                    Tamb_in     = [logFile.meteo.temperature];
                    Tamb_time   = [logFile.meteo.time];
                    [data.Thot(k), ~] = check_Thot_and_Tamb(Thot_time, Thot_in, Tamb_time, Tamb_in);
                    
                    

                    data.S_hot  = mean( [S_hot1(idx_hot1(k),:); S_hot2(idx_hot2(k),:)],1 ) ;
                    data.S_hot_diff  = S_hot1(idx_hot1(k),:) - S_hot2(idx_hot2(k),:) ;

                      
                    data.S_LN2 = mean( [S_LN2_1(idx_LN2_1(k),:); S_LN2_2(idx_LN2_2(k),:)],1 ) ;
                    data.S_LN2_diff = S_LN2_1(idx_LN2_1(k),:) - S_LN2_2(idx_LN2_2(k),:) ;


                    T_HOT_DIFF(k,:)=(data.S_hot_diff)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2);
                    T_LN2_DIFF(k,:)=(data.S_LN2_diff)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2);

                    % write data into struct
                    if nonempty_cold
                      data.S_cold = mean( [S_cold1(idx_cold1(k),:); S_cold2(idx_cold2(k),:)],1 ) ;
                      data.S_cold_diff = S_cold1(idx_cold1(k),:) - S_cold2(idx_cold2(k),:);
                      data.mirror_el_cold(k) = mirror_cold_d1(idx_cold1(k));
                      T_COLD(k,:)=(data.S_cold-data.S_LN2)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2)+T_LN2;
                      T_COLD_DIFF(k,:)=(data.S_cold_diff-data.S_LN2)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2)+T_LN2;

                    end
                    
                    if nonempty_line
                      data.S_line = mean( [S_line1(idx_line1(k),:); S_line2(idx_line2(k),:)],1 );
                      data.S_line_diff =  S_line1(idx_line1(k),:) - S_line2(idx_line2(k),:);
                      mirror_el_line(k,1) = mirror_line_d1(k);
                      T_LINE(k,:)=(data.S_line-data.S_LN2)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2)+T_LN2;
                      T_LINE_DIFF(k,:)=(data.S_line_diff)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2);
                    else
                      mirror_el_line(k,1) = nan;
                      T_LINE(k,:) = nan_rep;
                      T_LINE_DIFF(k,:) = nan_rep;
                    end

                    if nonempty_ref
                        data.S_ref  = mean( [S_ref1(idx_ref1(k),:); S_ref2(idx_ref2(k),:)],1 ) ;
                        data.S_ref_diff  = S_ref1(idx_ref1(k),:) - S_ref2(idx_ref2(k),:);
                        mirror_el_ref(k,1)  = mirror_ref_d1(idx_ref1(k));
                        T_REF(k,:)=(data.S_ref-data.S_LN2)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2)+T_LN2;
                        T_REF_DIFF(k,:)=(data.S_ref_diff)./(data.S_hot-data.S_LN2).*(data.Thot(k)-T_LN2);
                    else
                        T_REF(k,:) = nan_rep;
                        T_REF_DIFF(k,:)= nan_rep;
                        mirror_el_ref(k,1) = nan;
                    end

                    y=(data.S_hot)./(data.S_LN2);
                    T_rec(k,:) = (data.Thot(k)-y.*T_LN2)./(y-1);
                else
                    disp('Measurement mismatch: skipping this one')
                end
                a = 1;

            end
        else
            [T_LINE, T_COLD, T_REF, T_rec] = deal(nan_rep,nan_rep,nan_rep,nan_rep );
            [T_LINE_DIFF, T_COLD_DIFF, T_REF_DIFF, T_HOT_DIFF, T_LN2_DIFF] = deal(nan_rep,nan_rep,nan_rep,nan_rep, nan_rep );
            [mirror_el_line,mirror_el_ref, data.Thot,] = deal(nan,nan,nan) ;
        end

        calibratedSpectra(m_loop).TSys = mean(T_rec, 1 , 'omitnan');
        calibratedSpectra(m_loop).stdTSys = mean(std(T_rec, 'omitnan'), 'omitnan');
        calibratedSpectra(m_loop).TCold = mean(T_COLD,1, 'omitnan');
        calibratedSpectra(m_loop).THot  = mean(data.Thot, 'omitnan');
        calibratedSpectra(m_loop).TCold = mean(T_COLD,1, 'omitnan');
        calibratedSpectra(m_loop).TColdDiff = mean(T_COLD_DIFF,1, 'omitnan');
        calibratedSpectra(m_loop).TLN2Diff = mean(T_LN2_DIFF,1, 'omitnan');
        calibratedSpectra(m_loop).THotDiff = mean(T_HOT_DIFF,1, 'omitnan');

        calibratedSpectra(m_loop).Tb_line = mean(T_LINE,1, 'omitnan');
        calibratedSpectra(m_loop).stdTb_line = std(T_LINE,0,1, 'omitnan'); 
        calibratedSpectra(m_loop).Tb_ref = mean(T_REF,1, 'omitnan');
        calibratedSpectra(m_loop).Tb_ref_diff = mean(T_REF_DIFF,1, 'omitnan');

        calibratedSpectra(m_loop).meanAngleLine = mean(mirror_el_line,1, 'omitnan');
        calibratedSpectra(m_loop).meanAngleRef = mean(mirror_el_ref,1, 'omitnan');
  
        clear idx_ref idx_hot idx_cold
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