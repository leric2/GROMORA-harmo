function plot_integrated_spectra_generic(calibrationTool,integratedSpectra,lowerLim,upperLim)
% Just for a first look
try   
    % plotting a spectra every numberOfSpectraToGroup measurements
    TOD={};
    N=24;
    fig=figure('visible','off');
    clf
    set(gcf, 'PaperPosition', [1 1 19 27.7])
    orient(fig,'landscape')
    %suptitle([calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)])
    cm = colormap(parula(N));
    TOD{1}=[num2str(0) ' h'];
    hold on
    subplot(2,2,[1,3]); 
    for i=1:N
        plot(integratedSpectra(i).if,integratedSpectra(i).Tb,'Color',cm(i,:));
        %title([calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)])
        xlabel('IF [MHz]')
        ylabel('T_B [K]')
        
        ylim([lowerLim,upperLim])
        TOD{i}=num2str(integratedSpectra(i).TOD);
        hold on
    end
    title('Non Corrected, all channels')
    grid on
    legend(TOD,'Location','southoutside','NumColumns',4)
    
    subplot(2,2,2); 
    for i=1:N
        plot(integratedSpectra(i).if,integratedSpectra(i).TbWinCorr.*integratedSpectra(i).channelsQuality,'Color',cm(i,:));
        xlabel('IF [MHz]')
        ylabel('T_B [K]')
        
        ylim([lowerLim,upperLim])
        hold on
    end
    title('Corrected (windows), good channels')
    grid on
    
    subplot(2,2,4);
    for i=1:N
       
        plot(integratedSpectra(i).if,integratedSpectra(i).TbTroposphericWindowCorr.*integratedSpectra(i).channelsQuality,'Color',cm(i,:));
        
        xlabel('IF [MHz]')
        ylabel('T_B [K]')

        %ylim([lowerLim,upperLim])
        hold on
    end
    title('Corrected (troposphere), good channels')
    grid on
    
    
    % saveas(gcf,[retrievalTool.level1Folder 'calibratedHourlySpectra_' correctedSpectra.date],'jpg')
    print([calibrationTool.level1Folder calibrationTool.instrumentName '_integratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpdf','-fillpage')
    close

catch ME
    warning(ME.identifier,'Problem Plotting')
end

end

