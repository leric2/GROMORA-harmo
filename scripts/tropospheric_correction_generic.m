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
    
    Tmean = calibratedSpectra(t).meanAirTemperature - deltaT;
    
    % Ingold, v2: CHECK UNITS FOR T
    % T_mean  = (-18.772 + 0.7721 * calibratedSpectra(t).meanAirTemperature + 0.1452 * calibratedSpectra(t).meanRelHumidity);
    
    Tbg = 2.7;                   % Microwave background
    
    
    % Linear fit of spectrum's wings
    
    f_trop_corr  = [calibratedSpectra(t).freq(100:5000) calibratedSpectra(t).freq(end-5000:end-100)];
    
    % Clean brightness temperature (without spurious channels)
    Tb_temp=calibratedSpectra(t).Tb;
    
    Tb_temp(find(calibratedSpectra(t).channelsQuality==0))=NaN;
    
    % Not needed to interpolate and extrapolate
%     ib=find(isnan(Tb_temp));
%     ig=find(calibratedSpectra(t).channelsQuality==1);
%     

    
    Tb_trop_corr = [Tb_temp(100:5000) Tb_temp(end-5000:end-100)];
    ig=find(~isnan(Tb_trop_corr));
    ib=find(isnan(Tb_trop_corr));
    
    % If we find some nan in here, we interpolate it:
    if ~isempty(ig) & ~isempty(ib)
        Tb_trop_corr(ib)=interp1(f_trop_corr(ig),Tb_trop_corr(ig),f_trop_corr(ib),'linear');
    elseif isempty(ig)
        disp('Problem with this spectra : ');
        disp(calibratedSpectra(t).timeOfDay);
    end
        % Tb_clean(logical(calibratedSpectra(t).channelsQuality))=NaN;
    
    % Smooth the spectra 
    Tb_trop_corr_smoothed=conv(Tb_trop_corr,ones(11,1)/11,'same');
    [p,s,mu] = polyfit(f_trop_corr,Tb_trop_corr_smoothed,1 );  % linear fit
    Twing    = polyval(p, calibratedSpectra(t).freq, [], mu);        % polynomial evaluation
    
    
    % Transmitance calculated (Ingold)
    transmittance = (Tmean - Twing)./(Tmean - Tbg);
    
    % Troposph. corr.
    calibratedSpectra(t).TbTroposphericCorr = (Tb_temp - Tmean*(1-transmittance) ) ./ transmittance;
    calibratedSpectra(t).troposphericTransmittance = transmittance;
    calibratedSpectra(t).meanTroposphericTransmittance  = mean(transmittance);
    

end