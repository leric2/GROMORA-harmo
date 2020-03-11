function plot_spectra_generic(retrievalTool,drift,calibratedSpectra,lowerLim,upperLim,N)
% Just for a first look
try
% plotting a spectra every N measurements 
l=floor(linspace(1,length(calibratedSpectra),N));

%TOD={};

figure
clf
set(gcf, 'PaperPosition', [1 1 19 27.7])
suptitle([retrievalTool.dateStr(1:4) '-' retrievalTool.dateStr(6:7) '-' retrievalTool.dateStr(9:10)])
subplot(3,2,1); plot(drift.t, drift.Tn, 'k'), hold on, plot(drift.t,drift.TSysLog,'y'), ylabel('Tn [K]'), xlim([0,24]),ylim([nanmedian(drift.Tn)-80,nanmedian(drift.Tn)+80])
subplot(3,2,2); plot(drift.t, drift.Ta ,'g'), ylabel('Ta [K]'), xlim([0,24])
%set(gca, 'ColorOrder', [1 0.5 0.5; 0.2 0.2 0.2, 0 0 1],'NextPlot', 'replacechildren');
colors = {'r','g','b'};
subplot(3,2,3); 
for i=1:3
    plot(drift.t, drift.a(i,:),colors{i}), hold on, ylabel('Counts [-]'),xlim([0,24]),ylim([nanmedian(drift.a(3,:))-200,nanmedian(drift.a(1,:))+200])
end
subplot(3,2,4); plot(drift.t, drift.T, 'r'),  ylabel('T Hot  [K]'), ylim([mean(drift.T)-0.5 mean(drift.T)+0.5]),xlim([0,24])
subplot(3,2,5);
for i=1:N
    plot(calibratedSpectra(l(i)).if,calibratedSpectra(l(i)).Tb);
    
    %plot(calibratedSpectra(l(i)).Tb)
    %hold on
    
    %plot(calibratedSpectra(i).meanFromTbDownAll)
    %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
    %xlabel('f [GHz]')
    %xlim(1e-9*[calibratedSpectra(l(i)).freq(1),calibratedSpectra(l(i)).freq(end)])
    ylabel('T_B [K]')
    ylim([lowerLim,upperLim])
    %TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
    hold on
end

subplot(3,2,6);
for i=1:N
    plot(calibratedSpectra(l(i)).if,calibratedSpectra(l(i)).TN);
    ylabel('TN [K]')
    ylim([100,5000])
    hold on
end
for i=5:6 subplot(3,2,i); xlabel('IF [MHz]'); end
%legend(TOD)
print([retrievalTool.level1Folder 'calibratedSpectra_' retrievalTool.dateStr '_' retrievalTool.spectrometer],'-dpdf','-fillpage')
close
%saveas(gcf,[retrievalTool.level1Folder 'calibratedSpectra_' retrievalTool.dateStr '_' retrievalTool.spectrometer],'jpg')
%close

catch ME
    warning(ME.identifier,'Problem Plotting')
end

end

