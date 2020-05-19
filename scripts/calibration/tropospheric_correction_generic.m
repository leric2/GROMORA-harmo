function spectra = tropospheric_correction_generic(spectra,deltaT)
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
for t = 1:length(spectra)
    spectra(t).transmittanceMethod='Ingold_v1';
    if spectra(t).meanAirTemperature==-9999
        spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).if));
        spectra(t).troposphericTransmittance = -9999;
        spectra(t).troposphericOpacity=-9999;
        spectra(t).meanTroposphericTransmittance  = -9999;
        spectra(t).TbtropWinCorr = -9999;
        
    else
        % Mean tropospheric temperature
        Tmean = spectra(t).meanAirTemperature - deltaT;
        % Ingold, v2: CHECK UNITS FOR T
        % T_mean  = (-18.772 + 0.7721 * calibratedSpectra(t).meanAirTemperature + 0.1452 * calibratedSpectra(t).meanRelHumidity);
        
        Tbg = 2.7;                   % Microwave background
        
        % Linear fit of spectrum's wings
        f_trop_corr  = [spectra(t).freq(100:2000) spectra(t).freq(end-2000:end-100)];
        
        % Clean brightness temperature (without spurious channels)
        % Achtung, we use now the corrected window Tb
        if isfield(spectra,'TbCorr')
            Tb_temp=spectra(t).TbCorr;
        else 
            Tb_temp=spectra(t).Tb;
        end
        
        Tb_temp(find(spectra(t).channelsQuality==0))=NaN;
        
        % Not needed to interpolate and extrapolate
        %     ib=find(isnan(Tb_temp));
        %     ig=find(calibratedSpectra(t).channelsQuality==1);
        %
        Tb_trop_corr = [Tb_temp(100:2000) Tb_temp(end-2000:end-100)];
        ig=find(~isnan(Tb_trop_corr));
        ib=find(isnan(Tb_trop_corr));
        
        % If we find some nan in here, we interpolate it:
        if ~isempty(ig) & ~isempty(ib)
            Tb_trop_corr(ib)=interp1(f_trop_corr(ig),Tb_trop_corr(ig),f_trop_corr(ib),'linear');
        elseif isempty(ig)
            disp('Problem with this spectra : ');
            disp(spectra(t).TOD);
        end
        % Tb_clean(logical(calibratedSpectra(t).channelsQuality))=NaN;
        
        % Smooth the spectra
        Tb_trop_corr_smoothed=conv(Tb_trop_corr,ones(11,1)/11,'valid');
        %f_trop_corr=f_trop_corr(1:length(Tb_trop_corr_smoothed));
        
        %[p,s,mu] = polyfit(f_trop_corr,Tb_trop_corr_smoothed,1 );  % linear fit
        %Twing    = polyval(p, spectra(t).freq, [], mu);        % polynomial evaluation
        
        % Keeping the easy method for now
        Twing=nanmean(Tb_trop_corr_smoothed);
        
        % Transmitance calculated (Ingold)
        transmittance = (Tmean - Twing)./(Tmean - Tbg);
        
        % Method used in GROMOS and SOMORA
        r = read_datafile(sprintf('/home/esauvageat/Documents/GROSOM/RetrievalFiles/standard.spectrum.above.troposphere.aa', path) , 'matrix');
        T_strat = ( mean(r(1:2)) + mean(r((end-1):end)) ) /2;
        meanTWings=mean(Twing);
        
        %spectra(t).oldMeanTroposphericTransmittance=(meanTWings - Tmean)/(T_strat - Tmean);
        % Troposph. corr.
        spectra(t).TbTroposphericWindowCorr = (Tb_temp - Tmean*(1-transmittance) ) ./ transmittance;
        
        spectra(t).troposphericTransmittance = transmittance;
        spectra(t).troposphericOpacity=-log(transmittance);
        spectra(t).meanTroposphericTransmittance  = mean(transmittance);
    
end
end