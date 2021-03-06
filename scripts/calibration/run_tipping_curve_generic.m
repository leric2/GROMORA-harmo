function TC = run_tipping_curve_generic(rawSpectra, logFile, calibrationTool)
%==========================================================================
% NAME          | run_tipping_curve_generic.m
% TYPE          | function
% AUTHOR(S)     | Franziska Schranz, Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Performing tipping curve calibration. Depeding on the
%               | instrument, it can be needed or optionnal.
%               |
%               |
%               |
%               |
% ARGUMENTS     | INPUTS: 1. rawSpectra: Matrix of raw data (1 line per
%               |           cycle, channels as columns.
%               |         2. logFile: standardized log file
%               |         3. calibrationTool:
%               |               - instrumentName
%               |               - level1Folder, dateStr
%               |               - backgroundMWTb
%               |               - numberOfChannels
%               |               - TC
%               |               - indiceHot, indiceAntenna, indiceCold,
%               |                 indiceTC
%               |               - referenceTime
%               |               - tippingSize
%               |
%               |
%               | OUTPUTS: 1. TC: structure array containing all
%               |               tipping curve data for this day. It is
%               |               saved as sub-structure within logFile.
%               |
%               |
%               | CALLS:
%               |       GROMORA:
%               |           1. calibrationTool.get_tipping_curve_data(rawSpectra,logFile, calibrationTool)
%               |
%               |       MIWARA:
%               |           1. calibrationTool.get_tipping_curve_data(rawSpectra,logFile, calibrationTool)
%               |           2. find_tau_iteratively
%==========================================================================


%% get tipping curve data
try
    TC_data = calibrationTool.get_tipping_curve_data(rawSpectra,logFile, calibrationTool);
catch ME
    warning('no TC data found for this day');
    TC_data = struct();
end

if strcmp(calibrationTool.instrumentName,'MIAWARA-C')
    %% find opacity (tau)
    for k = 1:length(TC_data)
        
        % check if tipping calibration needs to be done for 2 polarisations
        if isfield(TC_data{1}  ,'s_tipping_pol1')
            [TC(1).tau(k), TC(1).Teff(k), TC(1).Trec_median(k), TC(1).quality(k), TC(1).offset(k), TC(1).niter(k), TC(1).time(k)] = find_tau_iteratively(TC_data{k},calibrationTool, TC_data{k}.s_tipping_pol1, TC_data{k}.s_hot_pol1, TC_data{k}.s_cold_pol1 );
            [TC(2).tau(k), TC(2).Teff(k), TC(2).Trec_median(k), TC(2).quality(k), TC(2).offset(k), TC(2).niter(k), TC(2).time(k)] = find_tau_iteratively(TC_data{k},calibrationTool, TC_data{k}.s_tipping_pol2, TC_data{k}.s_hot_pol2, TC_data{k}.s_cold_pol2 );
        else
            [TC.tau(k), TC.Teff(k), TC.Trec_median(k), TC.quality(k), TC.offset(k), TC.niter(k), TC.time(k)]  = find_tau_iteratively(TC_data{k},calibrationTool);
        end
    end
    %% write tau and Teff into a txt file
    
    if length(TC)==1
        T = table(TC.time', TC.tau',TC.Teff','VariableNames',{'time','tau','Teff'});
        writetable(T,[calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' calibrationTool.dateStr '_1']);
        
    else
        T = table(TC(1).time',TC(1).tau',TC(1).Teff','VariableNames',{'time','tau','Teff'});
        writetable(T,[calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' calibrationTool.dateStr '_1']);
        
        T = table(TC(2).time',TC(2).tau',TC(2).Teff','VariableNames',{'time','tau','Teff'});
        writetable(T,[calibrationTool.level1Folder calibrationTool.instrumentName '_tau_and_Teff_' calibrationTool.dateStr '_2']);
    end
    
else
    %% For instruments using hot-cold calibration scheme
    if ~isempty(fieldnames(TC_data)) && length(find(logFile.Tipping_Curve_active)) ~= length(find(logFile.Tipping_Curve_active & logFile.Position == calibrationTool.indiceTC))
        if isfield(logFile.meteo,'air_temperature')
            Teff = nanmean([logFile.meteo.air_temperature])-calibrationTool.TC.deltaT;
        else
            disp('we said, no meteo data found so lets make a guess for Tair (10 degC)');
            Teff = 283 - calibrationTool.TC.deltaT;
        end
        for i =1:length(TC_data)
            am = 1./sind(TC_data(i).skyAngle);
            tau = calibrationTool.TC.tauInitTC;
            it_max = calibrationTool.TC.maxIterTC;
            N_it = 0;
            off = 10;
            TC_data(i).converged = 1;
            while abs(off) > calibrationTool.TC.offsetTC
                % Computation of T_cold with this tau
                TCold=calibrationTool.backgroundMWTb.*exp(-tau)+Teff.*(1-exp(-tau));
                %disp(T_cold)
                % Computation of Tb_theta for each angle with this T_cold
                Tb_theta=TC_data(i).THot+(TC_data(i).THot-TCold).*(TC_data(i).sky-TC_data(i).hot)./(TC_data(i).hot-TC_data(i).cold);
                %disp(Tb_theta)
                
                % For every angle, compute the corresponding tau_theta
                tau_theta=log((Teff-calibrationTool.backgroundMWTb)./(Teff-Tb_theta));
                %disp(tau_theta)
                
                % Linear fit for tau_theta vs airmass factor
                fit=polyfit(am,tau_theta,1);
                
                % The new tau is the slope of the linear regression
                off = fit(2);
                tau = fit(1);
                N_it=N_it+1;
                
                if N_it>it_max
                    %tau=NaN;
                    %disp('Failed to converge')
                    TC_data(i).converged = 0;
                    break
                end
            end
            
            %TC_data(i).Tb_fromTCLoads = calibrationTool.TCold + (TC_data(i).THot - calibrationTool.TCold) .* (TC_data(i).sky - TC_data(i).cold)./(TC_data(i).hot - TC_data(i).cold);
            %tau_slant = log((Teff-calibrationTool.backgroundMWTb)./(Teff-TC_data(i).Tb_fromTCLoads));
            
            TC_data(i).tauIter = real(tau);
        end
    end
    TC = TC_data;
    
end