function plot_spectra_generic(calibrationTool, drift, meteoData,calibratedSpectra,lowerLim,upperLim,N)
% Just for a first look
% plotting a spectra every N measurements 
l=floor(linspace(1,length(calibratedSpectra),N));

%TOD={};

fig = figure('visible','off');
clf
set(gcf, 'PaperPosition', [1 1 19 27.7])
orient(fig,'landscape')

if ~isempty(drift)
    subplot(3,2,1); plot(drift.dateTime, drift.Tn, 'k'), hold on, plot(drift.dateTime,drift.TSysLog,'y'), ylabel('Tn [K]') ,ylim([nanmedian(drift.Tn)-80,nanmedian(drift.Tn)+80])
    subplot(3,2,2); plot(drift.dateTime, drift.Ta ,'g'), ylabel('Ta [K]')
    %set(gca, 'ColorOrder', [1 0.5 0.5; 0.2 0.2 0.2, 0 0 1],'NextPlot', 'replacechildren');
    colors = {'r','g','b'};
    subplot(3,2,3); 
    for i=1:3
        plot(drift.dateTime, drift.a(i,:),colors{i}), hold on, ylabel('Counts [-]'),ylim([nanmedian(drift.a(3,:))-200,nanmedian(drift.a(1,:))+200])
    end
    subplot(3,2,4); plot(drift.dateTime, drift.T, 'r'),  ylabel('T Hot  [K]'), ylim([mean(drift.T)-0.5 mean(drift.T)+0.5])
end
if ~isempty(meteoData)
    ax = subplot(3,2,5);
    yyaxis(ax(1),'left')
    plot(ax(1),[meteoData.dateTime], [meteoData.air_temperature] ,'r')
   % set(ax(1),'ylim', [0,300])
    set(ax(1),'xlim', [0,24])
    set(ax(1),'YColor','r');
    ylabel(ax(1),({'airT [K]'}))
    
    yyaxis(ax(1),'right')
    plot(ax(1),[meteoData.dateTime], [meteoData.air_pressure] ,'k')
    set(ax(1),'xlim', [0,24])
    set(ax(1),'YColor','k');
    ylabel(ax(1),({'airP [hPa]'}))
    
    ax2 = subplot(3,2,6);
    yyaxis(ax2(1),'left')
    plot(ax2(1), [meteoData.dateTime], [meteoData.precipitation] ,'k')
    %set(ax2(1),'ylim', [0,300])
    set(ax2(1),'xlim', [0,24])
    set(ax2(1),'YColor','k');
    ylabel(ax2(1),({'precipitation [?]'}))
    
    yyaxis(ax2(1),'right')
    plot(ax2(1),[meteoData.dateTime], [meteoData.rel_humidity] ,'b')
    set(ax2(1),'ylim', [0,100])
    set(ax2(1),'xlim', [0,24])
    set(ax2(1),'YColor','b');
    ylabel(ax2(1),({'relH [%]'}))
end
for i=1:6; subplot(3,2,i); grid on; end
print(fig,[calibrationTool.level1Folder calibrationTool.instrumentName '_calibratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpsc','-fillpage')
close

fig2 = figure('visible','off');
clf
set(gcf, 'PaperPosition', [.1 .1 0.5, 0.5])
orient(fig2,'landscape')
cm = colormap(parula(N));
subplot(1,2,1);
for i=1:N
    plot(calibratedSpectra(l(i)).if,calibratedSpectra(l(i)).Tb,'Color',cm(i,:));
    %plot(calibratedSpectra(l(i)).freq,calibratedSpectra(l(i)).Tb);
    
    %plot(calibratedSpectra(l(i)).Tb)
    TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
    
    %plot(calibratedSpectra(i).meanFromTbDownAll)
    %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
    %xlabel('f [GHz]')
    %xlim(1e-9*[calibratedSpectra(l(i)).freq(1),calibratedSpectra(l(i)).freq(end)])
    ylabel('T_B [K]')
    ylim([lowerLim,upperLim])
    %TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
    hold on
end

subplot(1,2,2);
for i=1:N
    plot(calibratedSpectra(l(i)).if,calibratedSpectra(l(i)).TN,'Color',cm(i,:));
    %plot(calibratedSpectra(l(i)).freq,calibratedSpectra(l(i)).T_rec);
    ylabel('TN [K]')
    ylim([100,5000])
    hold on
end
legend(TOD,'Location','southoutside','NumColumns',4)
for i=1:2; subplot(1,2,i); grid on, xlabel('IF [MHz]'); end


%legend(TOD)
%print([calibrationTool.level1Folder 'calibratedSpectra_' calibrationTool.dateStr '_' calibrationTool.spectrometer],'-dpdf','-fillpage')
print(fig2,[calibrationTool.level1Folder calibrationTool.instrumentName '_calibratedSpectra_' calibrationTool.spectrometer '_' calibrationTool.dateStr],'-dpsc','-append','-fillpage')
%saveas(gcf,[retrievalTool.level1Folder 'calibratedSpectra_' retrievalTool.dateStr '_' retrievalTool.spectrometer],'jpg')
close

