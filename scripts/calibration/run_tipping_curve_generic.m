function TC = run_tipping_curve_generic(rawSpectra, logFile, calibrationTool)


%% get tipping curve data

TC_data = calibrationTool.get_tipping_curve_data(rawSpectra,logFile, calibrationTool);



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

