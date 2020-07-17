function plot_integrated_spectra_generic(calibrationTool,integratedSpectra,lowerLim,upperLim)
% Just for a first look
try   
    % plotting a spectra every numberOfSpectraToGroup measurements
    %TOD={};
    N=length(integratedSpectra);
    fig=figure('visible','off');
    %fig=figure();
    clf
    set(gcf, 'PaperPosition', [1 1 19 27.7])
    orient(fig,'landscape')
    %suptitle([calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)])
    cm = colormap(parula(N));
    %TOD{1}=[num2str(0) ' h'];
    
    if ~(nanmedian([integratedSpectra.TSys])==-9999)
        yLowTN = nanmedian([integratedSpectra.TSys])-100;
        yUpTN = nanmedian([integratedSpectra.TSys])+100;
    else
        yLowTN = calibrationTool.TSysCenterTh - 1000;
        yUpTN = calibrationTool.TSysCenterTh + 1000;
    end
    
    ax = subplot(3,2,1);
    grid on
    yyaxis(ax(1),'left')
    plot([integratedSpectra.dateTime],[integratedSpectra.TSys],'k');
    hold on
    plot([integratedSpectra.dateTime], nanmedian([integratedSpectra.TSys])*[integratedSpectra.outlierCalib],'mx');
    set(ax(1),'ylim',[yLowTN,yUpTN])
    %set(ax(1),'xlim', [0,24])
    set(ax(1),'YColor','k');
    ylabel(ax(1),({'TN [K]'}))
    
    yyaxis(ax(1),'right')
    plot([integratedSpectra.dateTime],[integratedSpectra.THot],'r');
    set(ax(1),'ylim', [calibrationTool.THotTh-2,calibrationTool.THotTh+2])
    %set(ax(1),'xlim', [0,24])
    set(ax(1),'YColor','r');
    ylabel(ax(1),({'THot [K]'}))
    
    if ~(nanmedian([integratedSpectra.meanTb])==-9999)
        yLow = nanmedian([integratedSpectra.meanTb])-40;
        yUp = nanmedian([integratedSpectra.meanTb])+60;
    else
        yLow = 0;
        yUp = 280;
    end
    
    ax2 = subplot(3,2,2);
    yyaxis(ax2(1),'left')
    plot(ax2, [integratedSpectra.dateTime],[integratedSpectra.meanTb],'g');
    hold on
    plot(ax2, [integratedSpectra.dateTime], 0.9*yUp*[integratedSpectra.outlierCalib],'mx');
   % set(ax(1),'ylim', [0,300])
    set(ax2,'YColor','g');
    ylabel(ax2(1),({'Tb [K]'}))
    set(ax2,'ylim', [yLow,yUp])
    yyaxis(ax2(1),'right')
    plot([integratedSpectra.dateTime],[integratedSpectra.meanAirTemperature ],'r');
   % set(ax(1),'ylim', [0,300])
    %set(ax(1),'xlim', [0,24])
    set(ax2(1),'YColor','r');
    ylabel(ax2(1),({'airT [K]'}))
    grid on
%     yyaxis(ax2(1),'right')
%     plot([integratedSpectra.dateTime],[integratedSpectra.THot],'r');
%     %set(ax(1),'xlim', [0,24])
%     set(ax2(1),'YColor','r');
%     ylabel(ax2(1),({'THot [K]'}))
    
    %ax3 = subplot(3,2,3);

    
%     yyaxis(ax3(1),'right')
%     plot([integratedSpectra.dateTime],[integratedSpectra.meanAirPressure],'k');
%     set(ax(1),'xlim', [0,24])
%     set(ax3(1),'YColor','k');
%     ylabel(ax3(1),({'airP [hPa]'}))

    ax3 = subplot(3,2,3);
    yyaxis(ax3(1),'left')
    plot([integratedSpectra.dateTime],[integratedSpectra.numAntSpectra],'g');
    hold on
    plot([integratedSpectra.dateTime],[integratedSpectra.numColdSpectra],'b');
    plot([integratedSpectra.dateTime],[integratedSpectra.numHotSpectra],'r');
    %set(ax(1),'ylim', [0,300])
    %set(ax(1),'xlim', [0,24])
    set(ax3(1),'YColor','k');
    ylabel(ax3(1),({'# of spectra [-]'}))
    
    yyaxis(ax3(1),'right')
    plot([integratedSpectra.dateTime],[integratedSpectra.numberOfAveragedSpectra],'k');
    set(ax3(1),'ylim', [-0.1,6.1])
    set(ax3(1),'YColor','k');
    ylabel(ax3(1),({'# cal cycle [-]'}))
    grid on
    
    ax4 = subplot(3,2,4);
    yyaxis(ax4(1),'left')
    plot([integratedSpectra.dateTime],[integratedSpectra.rainAccumulation ],'k');
   % set(ax(1),'ylim', [0,300])
    %set(ax(1),'xlim', [0,24])
    set(ax4(1),'YColor','k');
    ylabel(ax4(1),({'precipitation [...]'}))
    
    yyaxis(ax4(1),'right')
    plot([integratedSpectra.dateTime],100*[integratedSpectra.meanRelativeHumidity],'b');
    set(ax4(1),'ylim', [0,100])
    set(ax4(1),'YColor','b');
    ylabel(ax4(1),({'RH [%]'}))
    grid on

%     for i=1:N
%         %plot(integratedSpectra(i).if,integratedSpectra(i).TbWinCorr.*integratedSpectra(i).channelsQuality,'Color',cm(i,:));
%         plot(integratedSpectra(i).dateTime,integratedSpectra(i).TSys);
%         %xlabel('IF [MHz]')
%         %ylabel('T_B [K]')
%         
%         %lim([lowerLim,upperLim])
%         hold on
%     end
%    title('Corrected (windows), good channels')

    subplot(3,2,5); 
    for i=1:N
        plot(integratedSpectra(i).if,integratedSpectra(i).Tb,'Color',cm(i,:));
        %title([calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)])
        xlabel('IF [MHz]')
        ylabel('T_B [K]')
        
        ylim([lowerLim,upperLim])
        %TOD{i}=num2str(integratedSpectra(i).TOD);
        hold on
    end
    title('Non Corrected, all channels')
    grid on
    %legend(TOD,'Location','southoutside','NumColumns',3)
    
    subplot(3,2,6);
    for i=1:N
       
        plot(integratedSpectra(i).if,integratedSpectra(i).TbTroposphericWindowCorr.*integratedSpectra(i).channelsQuality,'Color',cm(i,:));
        
        xlabel('IF [MHz]')
        ylabel('T_B [K]')

        %ylim([lowerLim,upperLim])
        
        hold on
    end
    grid on
    title('Corrected, good channels')
    
    % saveas(gcf,[retrievalTool.level1Folder 'calibratedHourlySpectra_' correctedSpectra.date],'jpg')
    print([calibrationTool.level1Folder calibrationTool.instrumentName '_integratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpdf','-fillpage')
    %print(fig,[calibrationTool.level1Folder calibrationTool.instrumentName '_integratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpsc','-fillpage')
    close

catch ME
    warning(ME.identifier,'%s',['Plotting integration problem: ' ME.message])
end

end

