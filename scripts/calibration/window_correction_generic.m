function spectra = window_correction_generic(calibrationTool,spectra)
%==========================================================================
% NAME          | 
% TYPE          |
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014
%               |
% ABSTRACT      |
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:
%               |
%               | OUTPUTS:
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================

% Window correction
for t = 1:length(spectra)
    freq=spectra(t).freq;
    if sum(spectra(t).Tb==-9999) == calibrationTool.numberOfChannels
        spectra(t).TbWinCorr = -9999*ones(1,calibrationTool.numberOfChannels);
    else
        if ~isnan(spectra(t).TWindow) 
            TWindow=spectra(t).TWindow;
        elseif ~isnan(spectra(t).meanAirTemperature)
            TWindow = spectra(t).meanAirTemperature;
        else
            error('no temperature found')
        end
            
        % Planck:
        TbWindowP= (2*calibrationTool.h*freq.^2)/(calibrationTool.lightSpeed^2)*(1)./(exp((t* freq)./(calibrationTool.kb*TWindow))-1);
    
        % Railey-Jeans ??
        % TbWindowRJ = (retrievalTool.lightSpeed^2 ./ (2*retrievalTool.kb*freq.^2) ) .* TbWindowP; 

        spectra(t).TbWinCorr  = (spectra(t).Tb - TbWindowP*(1-calibrationTool.tWindow))./calibrationTool.tWindow;
    end
end
end
