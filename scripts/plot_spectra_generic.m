function plot_spectra_generic(retrievalTool,calibratedSpectra,lowerLim,upperLim,N)
% Just for a first look

% plotting a spectra every N measurements 
l=floor(linspace(1,size(calibratedSpectra,2),N));

TOD={};

figure();
for i=1:N
    plot(calibratedSpectra(l(i)).Tb)
    %plot(calibratedSpectra(i).meanFromTbUpAll)
    %hold on
    %plot(calibratedSpectra(i).meanFromTbDownAll)
    %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
    ylim([lowerLim,upperLim])
    TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
    hold on
end
legend(TOD)

saveas(gcf,[retrievalTool.level1Folder 'calibratedSpectra_' calibratedSpectra(1).date],'jpg')
end

