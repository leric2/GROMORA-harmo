function spectra = tropospheric_correction_generic(spectra,calibrationTool,deltaT)
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
tropCorrType = calibrationTool.troposphericCorrection.type;

Tbg = calibrationTool.backgroundMWTb; % Microwave background

for t = 1:length(spectra)
    spectra(t).troposphericCorrType = tropCorrType;
    
    if spectra(t).mean_air_temperature==-9999 | sum(spectra(t).Tb == -9999) == calibrationTool.numberOfChannels
        spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
        spectra(t).troposphericTransmittance = -9999;
        spectra(t).troposphericOpacity=-9999;
        %spectra(t).meanTroposphericTransmittance  = -9999;
        %spectra(t).TbtropWinCorr = -9999;  
    else
        % Mean tropospheric temperature definition
        if strcmp(tropCorrType,'Ingold_v1') | strcmp(tropCorrType,'Ingold_v1_fit')
            Tmean = spectra(t).mean_air_temperature - deltaT;
        elseif strcmp(tropCorrType,'Ingold_v2')  
            %Ingold, v2: CHECK UNITS FOR T
            Tmean  = (-18.772 + 0.7721 * (spectra(t).mean_air_temperature-calibrationTool.zeroDegInKelvin) + 0.1452 * calibratedSpectra(t).mean_relative_humidity)+calibrationTool.zeroDegInKelvin;
        else
            error('Please specify a type of tropospheric correction to apply')
        end
        
        
        % Linear fit of spectrum's wings
        %f_trop_corr  = [spectra(t).freq(100:2000) spectra(t).freq(end-2000:end-100)];
        
        % Clean brightness temperature (without spurious channels)
        % Achtung, we use now the corrected WINDOW Tb
        if isfield(spectra,'TbWinCorr')
            Tb=spectra(t).TbWinCorr;
        else 
            Tb=spectra(t).Tb;
        end
        
        f_temp = spectra(t).frequencies;
        %Tb(isnan(spectra(t).channelsQuality))=NaN;
        
        Tb_temp = Tb;
        Tb_temp(isnan(spectra(t).channelsQuality)) = [];
        f_temp(isnan(spectra(t).channelsQuality)) = [];
        
        if length(Tb_temp) < 10*calibrationTool.troposphericCorrection.numberOfChannelsTropCorr
            spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
            spectra(t).troposphericTransmittance = -9999;
            spectra(t).troposphericOpacity=-9999;
            continue
        end 
        N = length(Tb_temp);
        skipChannels = calibrationTool.troposphericCorrection.skipFraction * N;
        
        lower = int16(skipChannels);
        upper = int16(skipChannels) + calibrationTool.troposphericCorrection.numberOfChannelsTropCorr;
        
        leftWingMeanFreq = nanmean(f_temp(lower:upper));
        leftWingMeanTb = nanmean(Tb_temp(lower:upper));
        rightWingMeanFreq = nanmean(f_temp(N-upper:N-lower));
        rightWingMeanTb = nanmean(Tb_temp(N-upper:N-lower));
        
        if strcmp(calibrationTool.troposphericCorrection.useWings,'both')
            Twing = (leftWingMeanTb+rightWingMeanTb)/2;
        elseif strcmp(calibrationTool.troposphericCorrection.useWings,'left')
            Twing = leftWingMeanTb;
        elseif strcmp(calibrationTool.troposphericCorrection.useWings,'right')
            Twing = rightWingMeanTb;
        else
            error('please specifiy the wings to use for tropospheric correction')
        end
        
        % Not needed to interpolate and extrapolate
        %     ib=find(isnan(Tb_temp));
        %     ig=find(calibratedSpectra(t).channelsQuality==1);
        %
%         Tb_trop_corr = [Tb_temp(100:2000) Tb_temp(end-2000:end-100)];
%         ig=find(~isnan(Tb_trop_corr));
%         ib=find(isnan(Tb_trop_corr));
%         
%         % If we find some nan in here, we interpolate it:
%         if ~isempty(ig) & ~isempty(ib)
%             Tb_trop_corr(ib)=interp1(f_trop_corr(ig),Tb_trop_corr(ig),f_trop_corr(ib),'linear');
%         elseif isempty(ig)
%             disp('Problem with this spectra : ');
%             disp(spectra(t).TOD);
%         end
        % Tb_clean(logical(calibratedSpectra(t).channelsQuality))=NaN;
        
        % Smooth the spectra
        %Tb_trop_corr_smoothed=conv(Tb_trop_corr,ones(11,1)/11,'valid');
        %f_trop_corr=f_trop_corr(1:length(Tb_trop_corr_smoothed));
        
        %[p,s,mu] = polyfit(f_trop_corr,Tb_trop_corr_smoothed,1 );  % linear fit
        %Twing    = polyval(p, spectra(t).freq, [], mu);        % polynomial evaluation
        
        % Keeping the easy method for now
        %Twing=nanmean(Tb_trop_corr_smoothed);
        
        
        if strcmp(tropCorrType,'Ingold_v1_fit')
            % frequency dependant opacity
            a = (rightWingMeanTb-leftWingMeanTb)/(rightWingMeanFreq - leftWingMeanFreq);
            b = rightWingMeanTb - a*rightWingMeanFreq;
            
            TbApprox = a * spectra(t).freq + b;
            transmittanceVector = (Tmean - TbApprox)./(Tmean - Tbg);
            
            if nanmean(transmittanceVector) > 0
                spectra(t).TbTroposphericWindowCorr = (Tb - Tmean*(1-transmittanceVector) ) ./ transmittanceVector;
                spectra(t).troposphericTransmittance = nanmean(transmittanceVector);
                spectra(t).troposphericOpacity=-log(nanmean(transmittanceVector));
            else
                spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
                spectra(t).troposphericTransmittance = -9999;
                spectra(t).troposphericOpacity=-9999;
            end
        else
            % Transmitance calculated (Ingold)
            transmittance = (Tmean - Twing)./(Tmean - Tbg);
        
            if transmittance > 0
                spectra(t).TbTroposphericWindowCorr = (Tb - Tmean*(1-transmittance) ) ./ transmittance;
                
                spectra(t).troposphericTransmittance = transmittance;
                spectra(t).troposphericOpacity=-log(transmittance);
                %spectra(t).meanTroposphericTransmittance  = mean(transmittance);
            else
                spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
                spectra(t).troposphericTransmittance = -9999;
                spectra(t).troposphericOpacity=-9999;
                %spectra(t).meanTroposphericTransmittance  = -9999;
                %spectra(t).TbtropWinCorr = -9999;
            end
        end
    end
end