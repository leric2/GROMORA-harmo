function plot_raw_spectra_generic(Spectra,lowerLim,upperLim,N)
%==========================================================================
% NAME          | plot_raw_spectra_generic.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Very quick plotting for raw data for the GROSOM project.
%               | 
%               | 
% ARGUMENTS     | INPUTS:   1. Spectra: raw spectra to plot
%               |           2. lowerLim
%               |           3. upperLim
%               |           4. N: number of raw spectra to plot
%               |
%               | OUTPUTS: -
%               |
%==========================================================================
try
    % plotting a spectra every N measurements
    l=floor(linspace(1,size(Spectra,1),N));
    
    figure();
    for i=1:N
        plot(Spectra(l(i),:))
        %plot(Spectra(10,:))
        ylim([lowerLim,upperLim])
        hold on
    end
catch ME
    warning(ME.identifier,'Problem with the plots of the raw spectra')
end

end
