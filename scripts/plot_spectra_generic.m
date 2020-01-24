function plot_spectra_generic(retrievalTool,calibratedSpectra,lowerLim,upperLim,N)
% Just for a first look
try
% plotting a spectra every N measurements 
l=floor(linspace(1,size(calibratedSpectra,2),N));

TOD={};

figure();
for i=1:N
    plot(calibratedSpectra(l(i)).freq/1e9,calibratedSpectra(l(i)).Tb)
    %plot(calibratedSpectra(i).meanFromTbUpAll)
    %hold on
    %plot(calibratedSpectra(i).meanFromTbDownAll)
    %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
    xlabel('f [GHz]')
    xlim(1e-9*[calibratedSpectra(l(i)).freq(1),calibratedSpectra(l(i)).freq(end)])
    ylabel('T_B [K]')
    ylim([lowerLim,upperLim])
    TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
    hold on
end
legend(TOD)
title([calibratedSpectra(1).date(1:4) '-' calibratedSpectra(1).date(6:7) '-' calibratedSpectra(1).date(9:10)])

saveas(gcf,[retrievalTool.level1Folder 'calibratedSpectra_' calibratedSpectra(1).date],'jpg')

% close

catch ME
    warning(ME.identifier,'Problem Plotting')
end

end

