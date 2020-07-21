function plot_spectra_generic(calibrationTool, drift, meteoData,calibratedSpectra,lowerLim,upperLim,N)
% Just for a first look
% plotting a spectra every N measurements 
try 
    l=floor(linspace(1,length(calibratedSpectra),N));
    
    TOD={};
    
    fig = figure('visible','off');
    %fig = figure();
    clf
    set(gcf, 'PaperPosition', [1 1 19 27.7])
    orient(fig,'landscape')
    
    xstart = calibrationTool.dateTime;
    xstop = calibrationTool.dateTime + day(1);
    
    limTNPlot = 2500;
    if ~isempty(drift)
        lowerLim = prctile(drift.Ta,2) - 20;
        upperLim = prctile(drift.Ta,98) + 20;
        
        limTNPlot = nanmedian(drift.Tn);
        
        ax = subplot(3,2,1); plot(ax,drift.dateTime, drift.Tn, 'k'), hold on, plot(ax,drift.dateTime,drift.TSysLog,'y'), ylabel('Tn [K]') , xlim([xstart,xstop]), ylim([limTNPlot-200,limTNPlot+200])
        plot(ax, [calibratedSpectra.meanAntTime], limTNPlot*[calibratedSpectra.outlierCalib]-100,'mx'), xlim([xstart,xstop]);
        
        ax2 = subplot(3,2,2); plot(ax2, drift.dateTime, drift.Ta ,'g'), hold on,ylabel('Ta [K]'), xlim([xstart,xstop]), ylim([0,nanmedian(drift.Ta)+120])
        plot(ax2, [calibratedSpectra.meanAntTime], nanmedian(drift.Ta)*[calibratedSpectra.outlierCalib],'mx'), xlim([xstart,xstop]);
        if ~isempty(drift.outlierCold) && length(drift.outlierCold) < 100
            plot(ax, drift.outlierCold,nanmedian(drift.Tn)-40,'bx'), xlim([xstart,xstop])
            plot(ax2, drift.outlierCold,10,'bx'), xlim([xstart,xstop])
        end
        if ~isempty(drift.outlierHot) && length(drift.outlierHot) < 100
            plot(ax, drift.outlierHot,nanmedian(drift.Tn)-50,'rx'), xlim([xstart,xstop])
            plot(ax2, drift.outlierHot,20,'rx'), xlim([xstart,xstop])
        end
        if ~isempty(drift.outlierSky) && length(drift.outlierSky) < 100
            plot(ax, drift.outlierSky,nanmedian(drift.Tn)-60,'gx'), xlim([xstart,xstop])
            plot(ax2, drift.outlierSky,30,'gx'), xlim([xstart,xstop])
        end
        %plot(logFile.dateTime(drift.outlierSky),200*ones(length(drift.outlierSky)),'x')
        %set(gca, 'ColorOrder', [1 0.5 0.5; 0.2 0.2 0.2, 0 0 1],'NextPlot', 'replacechildren');
        colors = {'r','g','b'};
        subplot(3,2,4);
        for i=1:3
            plot(drift.dateTime, drift.a(i,:),colors{i}), hold on, ylabel('Counts [-]'),ylim([prctile(drift.a(3,:),2)-200,prctile(drift.a(1,:),98)+200]), xlim([xstart,xstop])
        end
        %     if ~isempty(drift.outlierCold)
        %         for out = 1:length(drift.outlierCold)
        %             line([drift.outlierCold(out) drift.outlierCold(out)],[nanmedian(drift.a(3,:))-200,nanmedian(drift.a(1,:))+200],'Color','b');
        %         end
        %     end
        %     if ~isempty(drift.outlierSky)
        %         for out = 1:length(drift.outlierSky)
        %             line([drift.outlierSky(out) drift.outlierSky(out)],[nanmedian(drift.a(3,:))-200,nanmedian(drift.a(1,:))+200],'Color','g');
        %         end
        %     end
        %     if ~isempty(drift.outlierHot)
        %         for out = 1:length(drift.outlierHot)
        %             line([drift.outlierHot(out) drift.outlierHot(out)],[nanmedian(drift.a(3,:))-200,nanmedian(drift.a(1,:))+200],'Color','r');
        %         end
        %     end
        subplot(3,2,3); plot(drift.dateTime, drift.T, 'r'),  ylabel('T Hot  [K]'),xlim([xstart,xstop]), ylim([mean(drift.T)-0.5 mean(drift.T)+0.5])
    end
    if ~isempty(fieldnames(meteoData))
        ax = subplot(3,2,5);
        yyaxis(ax(1),'left')
        plot(ax(1),[meteoData.dateTime], [meteoData.air_temperature] ,'r')
        % set(ax(1),'ylim', [0,300])
        set(ax(1),'xlim',[xstart,xstop])
        set(ax(1),'YColor','r');
        ylabel(ax(1),({'airT [K]'}))
        
        yyaxis(ax(1),'right')
        plot(ax(1),[meteoData.dateTime], [meteoData.air_pressure] ,'k')
        set(ax(1),'xlim',[xstart,xstop])
        set(ax(1),'YColor','k');
        ylabel(ax(1),({'airP [hPa]'}))
        
        ax2 = subplot(3,2,6);
        yyaxis(ax2(1),'left')
        plot(ax2(1), [meteoData.dateTime], [meteoData.precipitation] ,'k')
        %set(ax2(1),'ylim', [0,300])
        set(ax2(1),'xlim',[xstart,xstop])
        set(ax2(1),'YColor','k');
        ylabel(ax2(1),({'precipitation [?]'}))
        
        yyaxis(ax2(1),'right')
        plot(ax2(1),[meteoData.dateTime], [meteoData.rel_humidity] ,'b')
        set(ax2(1),'ylim', [0,100])
        set(ax2(1),'xlim',[xstart,xstop])
        set(ax2(1),'YColor','b');
        ylabel(ax2(1),({'relH [%]'}))
    end
    
    %subplot(4,2,7); plot([calibratedSpectra.meanAntTime], [calibratedSpectra.flagged]);
    
    for i=1:6; subplot(3,2,i); grid on; end
    print(fig,[calibrationTool.level1Folder calibrationTool.instrumentName '_calibratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpsc','-fillpage')
    close
    
    fig2 = figure('visible','off');
    %fig2 = figure();
    clf
    set(gcf, 'PaperPosition', [.1 .1 0.5, 0.5])
    orient(fig2,'landscape')
    cm = colormap(parula(N));
    subplot(1,2,1);
    count=1;
    for i=1:N
        if ~(calibratedSpectra(l(i)).outlierCalib == 1)
            plot(calibratedSpectra(l(i)).if,calibratedSpectra(l(i)).Tb,'Color',cm(i,:));
            %plot(calibratedSpectra(l(i)).freq,calibratedSpectra(l(i)).Tb);
        
            %plot(calibratedSpectra(l(i)).Tb)
            
        
        %plot(calibratedSpectra(i).meanFromTbDownAll)
        %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
        %xlabel('f [GHz]')
        %xlim(1e-9*[calibratedSpectra(l(i)).freq(1),calibratedSpectra(l(i)).freq(end)])
            ylabel('T_B [K]')
        
            ylim([lowerLim,upperLim])
            TOD{count}=num2str(calibratedSpectra(l(i)).timeOfDay);
            count = count + 1;
        %TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
        end
        
        hold on
    end
    
    subplot(1,2,2);
    for i=1:N
        if ~(calibratedSpectra(l(i)).outlierCalib == 1)
            plot(calibratedSpectra(l(i)).if,calibratedSpectra(l(i)).TN,'Color',cm(i,:));
        %plot(calibratedSpectra(l(i)).freq,calibratedSpectra(l(i)).T_rec);
            ylabel('TN [K]')
            ylim([limTNPlot-1000,limTNPlot+1000])
        end
        hold on
    end
    
    
    legend(TOD,'Location','southoutside','NumColumns',4)
    for i=1:2; subplot(1,2,i); grid on, xlabel('IF [MHz]'); end
    
    
    %legend(TOD)
    %print([calibrationTool.level1Folder 'calibratedSpectra_' calibrationTool.dateStr '_' calibrationTool.spectrometer],'-dpdf','-fillpage')
    print(fig2,[calibrationTool.level1Folder calibrationTool.instrumentName '_calibratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpsc','-append','-fillpage')
    %saveas(gcf,[retrievalTool.level1Folder 'calibratedSpectra_' retrievalTool.dateStr '_' retrievalTool.spectrometer],'jpg')
    close
catch ME
    warning(ME.identifier,'%s',['Plotting calibration problem: ' ME.message])
end
