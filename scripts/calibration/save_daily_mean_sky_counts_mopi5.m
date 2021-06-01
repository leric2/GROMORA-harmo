function save_daily_mean_sky_counts_mopi5(calibrationTool, logFile, rawSpectra)

    skyInd = find(logFile.Position==calibrationTool.indiceAntenna & logFile.Measurement_NoiseDiode'==0); % Antenna
    coldInd = find(logFile.Position==calibrationTool.indiceCold & logFile.Measurement_NoiseDiode'==0); % Cold
    hotInd = find(logFile.Position==calibrationTool.indiceHot & logFile.Measurement_NoiseDiode'==0); % Hot

    daily_sky_counts = struct();
    daily_sky_counts.sky_counts = nanmean(rawSpectra(skyInd,:),1);
    daily_sky_counts.hot_counts = nanmean(rawSpectra(hotInd,:),1);
    daily_sky_counts.cold_counts = nanmean(rawSpectra(coldInd,:),1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency vector
    if calibrationTool.IQProcessing
        daily_sky_counts.if = calibrationTool.samplingRateFFTS/2 * [-1:2/calibrationTool.numberOfChannels:1-2/calibrationTool.numberOfChannels];
        
        daily_sky_counts.observationFreq=calibrationTool.observationFreq;
        
        daily_sky_counts.LOFreqTot=calibrationTool.LOFreqTot;
    
        daily_sky_counts.freq=daily_sky_counts.if*1e6+daily_sky_counts.LOFreqTot;
    
        daily_sky_counts.df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    else
        daily_sky_counts.if = calibrationTool.samplingRateFFTS/2 * [0:1/calibrationTool.numberOfChannels:1-1/calibrationTool.numberOfChannels];
    
        daily_sky_counts.observationFreq=calibrationTool.observationFreq;
        
        daily_sky_counts.LOFreqTot=calibrationTool.LOFreqTot;
    
        daily_sky_counts.freq=daily_sky_counts.if*1e6+daily_sky_counts.LOFreqTot;
    
        daily_sky_counts.df=calibrationTool.samplingRateFFTS/(2*calibrationTool.numberOfChannels);
    end
    
    save(['/home/es19m597/Documents/MOPI5/mean_daily_counts_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'daily_sky_counts')
end