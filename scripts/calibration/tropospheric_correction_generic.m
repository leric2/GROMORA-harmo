function spectra = tropospheric_correction_generic(spectra,calibrationTool)
%==========================================================================
% NAME          | tropospheric_correction_generic.m
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014, adapted 2020
%               |
% ABSTRACT      | Function doing tipping curve correction for GROSOM. It
%               | can be on a calibrated or integrated spectra and include
%               | different type of tropospheric correction.
%               |
%               |
% ARGUMENTS     | INPUTS:  1. spectra: standard spectra structure to
%               |            correct (calibrated or Integrated)
%               |          2. calibrationTool:
%               |            - troposphericCorrection: sub-structure
%               |            containing all needed parameters for a given
%               |            correction.
%               |            - backgroundMWTb
%               |            - deltaTCorr
%               |            - numberOfChannels
%               |            - zeroDegInKelvin
%               |            - referenceTime
%               |            - referenceTime
%               |
%               | OUTPUTS: 1. spectra with the following added fields:
%               |            - TbTroposphericWindowCorr
%               |            - troposphericTransmittance
%               |            - troposphericOpacity
%               |       
% COMMENTS      | The tropospheric corrected spectral Tb is done after the
%               | window correction and we therefore name it
%               | "TbTroposphericWindowCorr"
%               |
%==========================================================================

tropCorrType = calibrationTool.troposphericCorrection.type;

Tbg = calibrationTool.backgroundMWTb; % Microwave background

deltaT = calibrationTool.troposphericCorrection.deltaT;

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
            
            TbApprox = a * spectra(t).frequencies + b;
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
        if isfield(calibrationTool.logFile,'TC')
            TC = calibrationTool.logFile.TC;
            %isTC = isbetween([TC.dateTime],spectra(t).first_sky_time,spectra(t).theoreticalStartTime+minutes(calibrationTool.calibrationTime));
            isTC = ([TC.time] >= spectra(t).time_min & [TC.time] < spectra(t).last_sky_time);
            % check if there was a tc done during this cycle (only mean
            % datetime)
                
            % For now, we compute the TC tau only if there is 1 TC per
            % cycle. It means, that it computes tau only when the function
            % is applied on the calibratedSpectra and not on the integrated
            % Spectra. It can be changed...
            if sum(isTC) == 1
                TC(isTC).Tb = calibrationTool.TCold + (TC(isTC).THot_calib - calibrationTool.TCold) .* (TC(isTC).sky_spectra - TC(isTC).cold_spectra)./(TC(isTC).hot_spectra - TC(isTC).cold_spectra);
                TC(isTC).Tb_mean = nanmean(TC(isTC).Tb,2);
        
                tau_slant = log((Tmean-calibrationTool.backgroundMWTb)./(Tmean-TC(isTC).Tb_mean));
                am = 1./sind(TC(isTC).tipping_angle);
        
                % fit the airmass-slant opacity data pairs
                [p,s] = polyfit (am, tau_slant, 1);
                TC(isTC).tauCalibZenith = p(1);
        
                % In the beam direction:
                spectra(t).troposphericOpacityTC = TC(isTC).tauCalibZenith * 1/sind(spectra(t).mean_sky_elevation_angle);
            else
                spectra(t).troposphericOpacityTC = NaN;
            end
        else
            spectra(t).troposphericOpacityTC = NaN;
        end
    end
end