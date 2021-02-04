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
    else
        if ~isnan(spectra(t).TWindow) 
            TemperatureWindow = spectra(t).TWindow;
        else
            TemperatureWindow = calibrationTool.zeroDegInKelvin + 18;
            warning('no temperature found, correcting using standard window temperature (18°C)')
        end
        
        if calibrationTool.savePlanckIntensity
            % Planck radiation intensity
            BT_Planck = ((2*calibrationTool.h*frequencies.^3)/(calibrationTool.lightSpeed^2))./(exp((calibrationTool.h*frequencies)./(calibrationTool.kb*TemperatureWindow)) - 1);
            spectra(t).intensityPlanckWinCorr = (spectra(t).intensity_planck -  BT_Planck*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;
        end
        % Inverting attenuation from the window (see Fernandez, 2015)
        % might be working only using the RJ equivalent temperature ? TO
        % CHECK
        spectra(t).TbWinCorr  = (spectra(t).Tb -  TemperatureWindow*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;
        
        % How it was done before (wrong ?!):
        %spectra(t).TbWinCorr  = (spectra(t).Tb -  BT_Planck*(1-calibrationTool.transmittanceWindow))./calibrationTool.transmittanceWindow;

        % In practice, this would gives the exact same result:
        % spectra(t).TbWinCorr  = spectra(t).Tb./calibrationTool.transmittanceWindow;
        
        % It does not matter if window temperature is missing... unless the
        % above correction is not right....
    end
end
end
