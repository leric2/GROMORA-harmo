function spectra = window_correction_generic(calibrationTool,spectra)
%==========================================================================
% NAME          | window_correction_generic.m
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014, adapted 2020
%               |
% ABSTRACT      | Perform window correction for MWR
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:   1. calibrationTool:
%               |               - numberOfChannels
%               |               - h, lightSpeed, kb
%               |               - tWindow
%               |           2. spectra: a standard spectra (calibrated or
%               |               integrated) structure conatining at least a
%               |               set of frequencies and TB and ideally a
%               |               window temperature (otherwise, take the
%               |               air temperature.
%               |
%               | OUTPUTS:  1. spectra: including an extra field for the
%               |               corrected window corrected spectra
%               |               (TbWinCorr)
%               |
%==========================================================================

% Window correction
for t = 1:length(spectra)
    frequencies = spectra(t).frequencies;
    if sum(spectra(t).Tb==-9999) == calibrationTool.numberOfChannels
        spectra(t).TbWinCorr = -9999*ones(1,calibrationTool.numberOfChannels);
    else
        if ~isnan(spectra(t).TWindow) 
            TWindow = spectra(t).TWindow;
        else
            TWindow = calibrationTool.zeroDegInKelvin + 20;
            warning('no temperature found, correcting using standard window temperature (20°C)')
        end
        
        % Planck:
        TbWindowP = (2*calibrationTool.h*frequencies.^2)/(calibrationTool.lightSpeed^2)*(1)./(exp((t*frequencies)./(calibrationTool.kb*TWindow))-1);
    
        % Railey-Jeans ??
        % TbWindowRJ = (retrievalTool.lightSpeed^2 ./ (2*retrievalTool.kb*freq.^2) ) .* TbWindowP; 

        spectra(t).TbWinCorr  = (spectra(t).Tb - TbWindowP*(1-calibrationTool.tWindow))./calibrationTool.tWindow;
    end
end
end
