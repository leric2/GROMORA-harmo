function plot_spectra_generic(Spectra,lowerLim,upperLim,N)
% Just for a first look

% plotting a spectra every N measurements 
l=floor(linspace(1,size(Spectra,2),N));

figure();
for i=1:N
    plot(Spectra(l(i)).Tb)
    ylim([lowerLim,upperLim])
    hold on
end

end

