function plot_spectra_generic(calibratedSpectra,lowerLim,upperLim,N)
% Just for a first look

% plotting a spectra every N measurements 
l=floor(linspace(1,size(calibratedSpectra,2),N));

TOD={};

figure();
for i=1:N
    plot(calibratedSpectra(l(i)).Tb)
    %plot(calibratedSpectra(l(i)).Tb-calibratedSpectra2(l(i)).Tb)
    ylim([lowerLim,upperLim])
    TOD{i}=num2str(calibratedSpectra(l(i)).timeOfDay);
    hold on
end
legend(TOD)
end

