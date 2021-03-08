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
            warning('no brightness temperature corrected for the window found !')
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
        
        % For now, only works for both wings
        if calibrationTool.savePlanckIntensity
            if strcmp(calibrationTool.troposphericCorrection.useWings,'both')
                leftWingIntensity = nanmean(spectra(t).intensity_planck(lower:upper));
                rightWingIntensity = nanmean(spectra(t).intensity_planck(N-upper:N-lower));
                intensityWing = (leftWingIntensity+rightWingIntensity)/2;
                leftWingRJE = nanmean(spectra(t).TbRJEWinCorr(lower:upper));
                rightWingRJE = nanmean(spectra(t).TbRJEWinCorr(N-upper:N-lower));
                TwingRJE = (leftWingRJE+rightWingRJE)/2;
            else
                error('please, use both wings')
            end
        end
        
        if strcmp(calibrationTool.troposphericCorrection.useWings,'both')
            Twing = (leftWingMeanTb+rightWingMeanTb)/2;
            fwing = (leftWingMeanFreq+rightWingMeanFreq)/2;
        elseif strcmp(calibrationTool.troposphericCorrection.useWings,'left')
            Twing = leftWingMeanTb;
            fwing = leftWingMeanFreq;
        elseif strcmp(calibrationTool.troposphericCorrection.useWings,'right')
            Twing = rightWingMeanTb;
            fwing = rightWingMeanFreq;
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
            
            TbgCorr = Tbg + calibrationTool.h*spectra(t).frequencies/(2*calibrationTool.kb);
            %TmeanCorr = Tmean + calibrationTool.h*spectra(t).frequencies/(2*calibrationTool.kb);
            transmittanceVectorCorr = (Tmean - TbApprox)./(Tmean - TbgCorr);
            if nanmean(transmittanceVector) > 0
                spectra(t).TbTroposphericWindowCorr = (Tb - Tmean*(1-transmittanceVector) ) ./ transmittanceVector;
                spectra(t).troposphericTransmittance = nanmean(transmittanceVector);
                spectra(t).troposphericOpacity=-log(nanmean(transmittanceVector));
                spectra(t).troposphericOpacityFromPhysicalTemp=-log(nanmean(transmittanceVector));
                spectra(t).troposphericOpacityFromPhysicalTempCorr=-log(nanmean(transmittanceVectorCorr));
            else
                spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
                spectra(t).troposphericTransmittance = -9999;
                spectra(t).troposphericOpacity=-9999;
                spectra(t).troposphericOpacityFromPhysicalTemp=-9999;
                spectra(t).troposphericOpacityFromPhysicalTempCorr=-9999;
            end
            
            if calibrationTool.savePlanckIntensity
                % frequency dependant opacity with RJE brightness
                % temperature
                a = (rightWingRJE-leftWingRJE)/(rightWingMeanFreq - leftWingMeanFreq);
                b = rightWingRJE - a*rightWingMeanFreq;
                
                TmeanRJvector = rayleigh_jeans_equivalent_Tb(calibrationTool, Tmean, spectra(t).frequencies);
                TbgRJvector = rayleigh_jeans_equivalent_Tb(calibrationTool, Tbg, spectra(t).frequencies);
                TbRJEApprox = a * spectra(t).frequencies + b;
                
                % we get the same transmittance if we use intensities.
                transmittanceVectorRJE = (TmeanRJvector - TbRJEApprox)./(TmeanRJvector - TbgRJvector);
                
                if nanmean(transmittanceVectorRJE) > 0
                    spectra(t).intensityPlanckTropWinCorr = (spectra(t).intensityPlanckWinCorr - planck_function(calibrationTool, Tmean, spectra(t).frequencies).*(1-transmittanceVectorRJE) ) ./ transmittanceVectorRJE;
                    spectra(t).TbRJETropWinCorr = (spectra(t).TbRJEWinCorr - TmeanRJvector.*(1-transmittanceVectorRJE) ) ./ transmittanceVectorRJE;
                    
                    % this is the variable that we keep ?!!
                    spectra(t).TbTroposphericWindowCorr = (Tb - Tmean.*(1-transmittanceVectorRJE) ) ./ transmittanceVectorRJE;
                    
                    %spectra(t).TbRJETroposphericWinCorr2 = spectra(t).intensityPlanckTropWinCorr*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*(spectra(t).frequencies.^2));
                    
                    spectra(t).troposphericTransmittance = nanmean(transmittanceVectorRJE);
                    spectra(t).troposphericOpacity = -log(nanmean(transmittanceVectorRJE));
                else
                    spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
                    spectra(t).troposphericTransmittance = -9999;
                    spectra(t).troposphericOpacity=-9999;
                    
                end
            end
        else
            % Transmitance calculated (Ingold) with Planck temperature (-->
            % 2.736K for background. This is an Approximation !
            transmittance = (Tmean - Twing)./(Tmean - Tbg);
            %TbgCorr = Tbg + calibrationTool.h*spectra(t).frequencies/(2*calibrationTool.kb);
            %TmeanCorr = Tmean + calibrationTool.h*spectra(t).frequencies/(2*calibrationTool.kb);
            %transmittanceVectorCorr = (TmeanCorr - TbApprox)./(TmeanCorr - TbgCorr);
            if transmittance > 0
                spectra(t).TbTroposphericWindowCorr = (Tb - Tmean*(1-transmittance) ) ./ transmittance;
                spectra(t).troposphericTransmittance = transmittance;
                spectra(t).troposphericOpacity=-log(transmittance);
                spectra(t).troposphericOpacityFromPhysicalTemp=-log(transmittance);
                spectra(t).troposphericOpacityFromPhysicalTempCorr=-9999;
                %spectra(t).troposphericOpacityIntensity=-log(transmittanceIntensity);
                %spectra(t).meanTroposphericTransmittance  = mean(transmittance);
            else
                
                spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
                spectra(t).troposphericTransmittance = -9999;
                spectra(t).troposphericOpacity=-9999;
                spectra(t).troposphericOpacityFromPhysicalTemp=-9999;
                spectra(t).troposphericOpacityFromPhysicalTempCorr=-9999;
            end
            % More accurate value with RJE OR Intensity calibrated value
            if calibrationTool.savePlanckIntensity
                % The right way of computing the transmittance
                %transmittanceIntensity = (B_PlanckMean - intensityWing)./(B_PlanckMean - planck_function(calibrationTool, Tbg, fwing));
                % it is equal to:
                TmeanRJE = rayleigh_jeans_equivalent_Tb(calibrationTool, Tmean, fwing);
                TbgRJE = rayleigh_jeans_equivalent_Tb(calibrationTool, Tbg, fwing);
                
                % here, we only compute a single transmittance value
                transmittanceRJE = (TmeanRJE - TwingRJE)./(TmeanRJE - TbgRJE);
                
                if transmittanceRJE > 0
                    % to check if that makes sense
                    spectra(t).intensityPlanckTropWinCorr = (spectra(t).intensityPlanckWinCorr - planck_function(calibrationTool, Tmean, fwing)*(1-transmittanceRJE) ) ./ transmittanceRJE;
                    %spectra(t).TbPlanckTropWinCorr = (calibrationTool.h*spectra(t).frequencies/calibrationTool.kb)./log((2*calibrationTool.h*spectra(t).frequencies.^3)./(spectra(t).intensityPlanckTropWindowCorr *calibrationTool.lightSpeed^2) + 1);
                    spectra(t).TbRJETropWinCorr = (spectra(t).TbRJEWinCorr - TmeanRJE*(1-transmittanceRJE) ) ./ transmittanceRJE;
                    
                    spectra(t).TbTroposphericWindowCorr = (Tb - Tmean*(1-transmittanceRJE) ) ./ transmittanceRJE;
                    %I_RJE =  spectra(t).TbRJETroposphericWindcorr.*(2*calibrationTool.kb*spectra(t).frequencies.^2)/calibrationTool.lightSpeed^2;
                    %spectra(t).TbPlanckTropWinCorr2 = (calibrationTool.h*spectra(t).frequencies./calibrationTool.kb)./log((2*calibrationTool.h*spectra(t).frequencies.^3)./(I_RJE*calibrationTool.lightSpeed^2) + 1);
                    
                    spectra(t).troposphericTransmittance = transmittanceRJE;
                    spectra(t).troposphericOpacity = -log(transmittanceRJE);
                else
                    spectra(t).TbTroposphericWindowCorr = -9999*ones(1,length(spectra(1).intermediate_freq));
                    spectra(t).troposphericTransmittance = -9999;
                    spectra(t).troposphericOpacity=-9999;
                end
                
                
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
            % Also, we do not account for the window yet !
            if sum(isTC) == 1
                % Using physical temperatures
                TC(isTC).Tb = calibrationTool.TCold + (TC(isTC).THot_calib - calibrationTool.TCold) .* (TC(isTC).sky_spectra - TC(isTC).cold_spectra)./(TC(isTC).hot_spectra - TC(isTC).cold_spectra);
                
                % In terms of RJE:
                TbColdRJ = rayleigh_jeans_equivalent_Tb(calibrationTool, calibrationTool.TCold, TC(isTC).frequency);
                TbHotRJ = rayleigh_jeans_equivalent_Tb(calibrationTool, TC(isTC).THot_calib, TC(isTC).frequency);
                
                TC(isTC).TbRJE = TbColdRJ + (TbHotRJ - TbColdRJ) .* (TC(isTC).sky_spectra - TC(isTC).cold_spectra)./(TC(isTC).hot_spectra - TC(isTC).cold_spectra);
                
                TC(isTC).TbRJE_mean = nanmean(TC(isTC).TbRJE,2);
                TC(isTC).Tb_mean = nanmean(TC(isTC).Tb,2);
                
                TmeanRJE = planck_function(calibrationTool, Tmean, fwing)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*(fwing^2));
                TbgRJE = planck_function(calibrationTool, Tbg, fwing)*(calibrationTool.lightSpeed^2)./(2*calibrationTool.kb*(fwing^2)); 
                
                tau_slantRJ = log((TmeanRJE-TbgRJE)./(TmeanRJE-TC(isTC).TbRJE_mean));
                tau_slant = log((Tmean-calibrationTool.backgroundMWTb)./(Tmean-TC(isTC).Tb_mean));
                am = 1./sind(TC(isTC).tipping_angle);
                
                % fit the airmass-slant opacity data pairs TbRJE
                [pRJE,sRJE] = polyfit (am, tau_slantRJ, 1);
                TC(isTC).tauCalibZenithRJE = pRJE(1);
                
                % fit the airmass-slant opacity data pairs 
                [p,s] = polyfit (am, tau_slant, 1);
                TC(isTC).tauCalibZenith = p(1);
                
                % Difference between the 2 is minimal...
                spectra(t).tauZenithTC = TC(isTC).tauCalibZenith;
                spectra(t).tauZenithTCRJE = TC(isTC).tauCalibZenithRJE;
        
                % In the beam direction:
                spectra(t).troposphericOpacityTC = TC(isTC).tauCalibZenithRJE * 1/sind(spectra(t).mean_sky_elevation_angle);               
            else
                spectra(t).troposphericOpacityTC = NaN;
            end
        else
            spectra(t).troposphericOpacityTC = NaN;
        end
    end
end