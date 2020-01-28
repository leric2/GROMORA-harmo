function calibratedSpectra = tropospheric_correction_generic(calibratedSpectra,deltaT)
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

for t = 1:length(calibratedSpectra)

    % Mean tropospheric temperature
    
    Tmean = calibratedSpectra(t).meanAirTemperature+273.15 - deltaT;
    
    Tbg = 2.7;                   % Microwave background
    
    
    % Linear fit of spectrum's wings
    
    f_trop_corr  = [calibratedSpectra(t).freq(100:5000) calibratedSpectra(t).freq(end-5000:end-100)];
    
    % Clean brightness temperature (without spurious channels)
    Tb_clean=calibratedSpectra(t).Tb;
    % Tb_clean(logical(calibratedSpectra(t).channelsQuality))=NaN;
    
    Tb_trop_corr = [Tb_clean(100:5000) Tb_clean(end-5000:end-100)];
    
    [p,s,mu] = polyfit(f_trop_corr,Tb_trop_corr,1 );  % linear fit
    Twing    = polyval(p, calibratedSpectra(t).freq, [], mu);        % polynomial evaluation
    
    
    % Transmitance calculated (Ingold)
    
    transmittance = (Tmean - Twing)./(Tmean - Tbg);
    
    % Troposph. corr.
    
    calibratedSpectra(t).TbTroposphericCorr = ( calibratedSpectra(t).Tb - Tmean*(1-transmittance) ) ./ transmittance;
    calibratedSpectra(t).troposphericTransmittance = transmittance;
    calibratedSpectra(t).meanTroposphericTransmittance  = mean(transmittance);

end