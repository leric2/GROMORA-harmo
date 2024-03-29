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
%               |               - transmittanceWindow
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
        spectra(t).intensityPlanckWinCorr = -9999*ones(1,calibrationTool.numberOfChannels);
        spectra(t).Tb_RJE = -9999*ones(1,calibrationTool.numberOfChannels);
        spectra(t).TbRJEWinCorr = -9999*ones(1,calibrationTool.numberOfChannels);
        spectra(t).TbWinCorrPhysicalTemperature = -9999*ones(1,calibrationTool.numberOfChannels);
    else
        if ~isnan(spectra(t).TWindow) 
            TemperatureWindow = spectra(t).TWindow;
        elseif ~isnan(spectra(t).TRoom)
            TemperatureWindow = spectra(t).TRoom;
            if t==1
                warning('no temperature found, correcting using Room temperature')
            end
        else
            TemperatureWindow = calibrationTool.zeroDegInKelvin + 20;
            if t==1
                warning('no temperature found, correcting using standard window temperature (20°C)')
            end
        end
        
        if calibrationTool.savePlanckIntensity
            % Correcting with intensity:
            intensityWindow = planck_function(calibrationTool, TemperatureWindow, frequencies);
            spectra(t).intensityPlanckWinCorr = (spectra(t).intensity_planck -  intensityWindow*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;
            spectra(t).intensityPlanckWinCorr(spectra(t).intensityPlanckWinCorr < 0) = nan;
            spectra(t).TbWinCorr = planck_Tb(calibrationTool, spectra(t).intensityPlanckWinCorr, frequencies);
            
            % Correcting with RJE Tb (should be the same)
            TemperatureWindowRJE = rayleigh_jeans_equivalent_Tb(calibrationTool, TemperatureWindow, frequencies);
            spectra(t).Tb_RJE = rayleigh_jeans_equivalent_Tb(calibrationTool, spectra(t).Tb, frequencies);
            spectra(t).TbRJEWinCorr = (spectra(t).Tb_RJE -  TemperatureWindowRJE*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;
            %spectra(t).intensityFromRJE = spectra(t).TbRJEWinCorr.*(2*calibrationTool.kb.*frequencies.^2)./calibrationTool.lightSpeed^2;
            
            % Correcting with physical temperature (not accurate)
            spectra(t).TbWinCorrPhysicalTemperature  = (spectra(t).Tb -  TemperatureWindow*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;
            
            % works and gives the same...
            %I_RJE =  spectra(t).TbRJEWinCorr .*(2*calibrationTool.kb*frequencies.^2)/calibrationTool.lightSpeed^2;
            %spectra(t).TbPlanck2 = (calibrationTool.h*frequencies/calibrationTool.kb)./log((2*calibrationTool.h*frequencies.^3)./(I_RJE*calibrationTool.lightSpeed^2) + 1);
        else
            % Inverting attenuation from the window (see Fernandez, 2015)
            % might be working only using the RJ equivalent temperature ? TO
            % CHECK
            spectra(t).TbWinCorr  = (spectra(t).Tb -  TemperatureWindow*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;
        end
        % How it was done before (wrong ?!):
        %spectra(t).TbWinCorr  = (spectra(t).Tb -  BT_Planck*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;

        % In practice, this would gives the exact same result:
        % spectra(t).TbWinCorr  = spectra(t).Tb./calibrationTool.transmittanceWindow;
        
        % It does not matter if window temperature is missing... unless the
        % above correction is not right....
    end
end
end
