function plot_raw_spectra_generic(rawSpectra,upperLim,N)
% Just for a first look

% plotting a spectra every N measurements 
l=floor(linspace(1,size(rawSpectra,1),N));

figure();
for i=1:N
    plot(rawSpectra(l(i),:))
    ylim([0,upperLim])
    hold on
end

end

