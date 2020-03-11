function plot_hourly_spectra_generic(retrievalTool,integratedSpectra,lowerLim,upperLim)
% Just for a first look
try   
    % plotting a spectra every numberOfSpectraToGroup measurements
    TOD={};
    N=24;
    fig=figure();
    clf
    set(gcf, 'PaperPosition', [1 1 19 27.7])
    suptitle([retrievalTool.dateStr(1:4) '-' retrievalTool.dateStr(6:7) '-' retrievalTool.dateStr(9:10)])
    
    TOD{1}=[num2str(0) ' h'];
    hold on
    subplot(1,2,1); 
    for i=1:N
        plot(integratedSpectra(i).if,integratedSpectra(i).Tb);
        title('')
        xlabel('IF [MHz]')
        ylabel('T_B [K]')
        
        ylim([lowerLim,upperLim])
        TOD{i}=num2str(integratedSpectra(i).TOD);
        hold on
    end
    title('Non Corrected, all channels')
    grid on
    subplot(1,2,2); 
    
    for i=1:N
        plot(integratedSpectra(i).if,integratedSpectra(i).Tbcorr.*integratedSpectra(i).channelsQuality);
        xlabel('IF [MHz]')
        ylabel('T_B [K]')
        
        ylim([lowerLim,upperLim])
        hold on
    end
    title('Corrected, good channels')
    %legend(TOD,'Location','southoutside')
    grid on
    orient(fig,'landscape')
    % saveas(gcf,[retrievalTool.level1Folder 'calibratedHourlySpectra_' correctedSpectra.date],'jpg')
    print([retrievalTool.level1Folder 'integratedSpectra' retrievalTool.dateStr '_' retrievalTool.spectrometer],'-dpdf','-fillpage')
    close

catch ME
    warning(ME.identifier,'Problem Plotting')
end

end

