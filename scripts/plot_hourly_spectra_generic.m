function plot_hourly_spectra_generic(retrievalTool,correctedSpectra,lowerLim,upperLim)
% Just for a first look
try
    numberOfSpectraToGroup=60/retrievalTool.calibrationTime;
    
    N=24;
    
    Tb=ones(length(correctedSpectra),retrievalTool.numberOfChannels)*NaN;
    for t = 1:length(correctedSpectra)
        Tb(t,:)=correctedSpectra(t).Tb;
    end
    
    % plotting a spectra every numberOfSpectraToGroup measurements
    TOD={};
    
    figure();
    
    meanHourlySpectra = mean(Tb(1:6,:),1);
    plot(correctedSpectra(1).freq/1e9,meanHourlySpectra)
    TOD{1}=[num2str(0) ' h'];
    hold on
    for i=2:N
        meanHourlySpectra= mean(Tb((i-1)*numberOfSpectraToGroup:i*numberOfSpectraToGroup,:),1);
        
        plot(correctedSpectra(1).freq/1e9,meanHourlySpectra)
        %plot(calibratedSpectra(i).meanFromTbUpAll)
        %hold on
        %plot(calibratedSpectra(i).meanFromTbDownAll)
        %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
        xlabel('f [GHz]')
        %xlim(1e-9*[correctedSpectra(1).freq(1),correctedSpectra.freq(end)])
        ylabel('T_B [K]')
        ylim([lowerLim,upperLim])
        TOD{i}=[num2str(i-1) ' [h]'];
        hold on
    end
    legend(TOD)
    
    saveas(gcf,[retrievalTool.level1Folder 'calibratedHourlySpectra_' correctedSpectra.date],'jpg')

% close

catch ME
    warning(ME.identifier,'Problem Plotting')
end

end

