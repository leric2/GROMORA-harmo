function [calibrationTool, level1] = run_integration(calibrationTool)
%==========================================================================
% NAME      | run_integration.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 01.2020
%           |
% ABSTRACT  | The main function executing the integration part for the
%           | instrument defined in calibrationTool. Some parts and
%           | functions are dependent on the instrument that we want to
%           | integrate and are stored in calibrationTool. 
%           | 
%           |
% ARGUMENTS | INPUTS: 1. calibrationTool: structure containing all
%           | information about the integration we want to perform.
%           | Documentation about this structure can be found in external
%           | document.
%           |
%           | OUTPUTS: 1. calibrationTool
%           |          2. level1: structure containing both level 1a and 1b
%           |           
% SAVE      | level1b netCDF file and plots
%           |
% CALLS     | Some depends on instruments, all are stored in calibrationTool:
%           | read_level1a
%           | add_meteo_data
%           | check_channel_quality (2x)
%           | window_correction (2x)
%           | tropospheric_correction_generic (2x)
%           | integrate_calibrated_spectra
%           | check_integrated
%           | plot_integrated_spectra
%           | save_level1b
%           |
% COMMENTS  | External documentation can be found for this function on the
%           | git server of IAP.
%==========================================================================
% Check that the calibration file (level1a) exist for this day
%     if ~isfield(calibrationTool, 'filenameLevel1a')
%         calibrationTool.filenameLevel1a=[calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr calibrationTool.extraName '.nc'];
%     end
% if ~exist(calibrationTool.filenameLevel1a,'file')
%     error('No calibration data found for this day')
%  end

% Start integration
disp(['Starting the integration process for ' calibrationTool.instrumentName ' ' calibrationTool.spectrometer ': ' calibrationTool.dateStr])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading the level1a data and adding the meteo data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We create level1b structure which will contain both the calibrated
% spectra read from the level1a data (level1b.calibratedSpectra and the 
% integrated spectra that will be added later (level1b.integration) 
level1 = struct();

[level1.calibratedSpectra, meteoData, calibrationTool] = calibrationTool.read_level1a(calibrationTool);

if isnan(meteoData.dateNum)
    %TODO
    disp('no meteo data for this day...')
    disp('Deactivate tropospheric transmittance filtering')
    calibrationTool.missing_meteo = true;
    level1.calibratedSpectra = calibrationTool.add_meteo_data(calibrationTool, meteoData, level1.calibratedSpectra);
else
    calibrationTool.missing_meteo = false;
    level1.calibratedSpectra = calibrationTool.add_meteo_data(calibrationTool, meteoData, level1.calibratedSpectra);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction and integration of the calibrated spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the quality of the channels and flagging the potential bad ones
% (we do not remove any)
filterType = calibrationTool.filterTypeChannelQualityCal;
level1.calibratedSpectra = calibrationTool.check_channel_quality(level1.calibratedSpectra,calibrationTool,filterType);

% Performing window correction
level1.calibratedSpectra = calibrationTool.window_correction(calibrationTool,level1.calibratedSpectra);

% Compute tropospheric transmittance and correction for every calibrated
% spectra.
level1.calibratedSpectra = calibrationTool.tropospheric_correction(level1.calibratedSpectra, calibrationTool);

% Integrating the "good spectra" based on tropospheric transmittance and
% calibration flags. --> To improve. Maybe introduce weighted mean of
% spectra based on tropospheric transmittance ?
level1.integratedSpectra = calibrationTool.integrate_calibrated_spectra(calibrationTool,level1.calibratedSpectra);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction and check of the integrated spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now on the integrated spectra; check the quality of the channels and 
% flagging the potential bad ones (we do not remove any).
filterType = calibrationTool.filterTypeChannelQualityInt;
level1.integratedSpectra = calibrationTool.check_channel_quality(level1.integratedSpectra,calibrationTool,filterType);

% Performing window correction
level1.integratedSpectra = calibrationTool.window_correction(calibrationTool,level1.integratedSpectra);

% Compute tropospheric transmittance and correction for every integrated
% spectra.
level1.integratedSpectra = calibrationTool.tropospheric_correction(level1.integratedSpectra, calibrationTool);

% Check integrated spectra and define the flags for level1b
level1.integratedSpectra = calibrationTool.check_integrated(calibrationTool, level1.integratedSpectra);

% sideband correction ?
% TODO

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting and saving level 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting integrated and corrected spectra
if calibrationTool.integratedSpectraPlot && nansum([level1.integratedSpectra.outlierCalib])<length(level1.integratedSpectra) && ~calibrationTool.missing_meteo
    calibrationTool.plot_integrated_spectra(calibrationTool,level1.integratedSpectra)
elseif calibrationTool.integratedSpectraPlot && length(level1.integratedSpectra) > 0
    calibrationTool.plot_integrated_spectra(calibrationTool,level1.integratedSpectra)
end

if calibrationTool.saveLevel1b
    % Saving integrated spectra (level1b) into NetCDF-4 file
    disp('Saving Level 1b...')
    calibrationTool  =  calibrationTool.save_level1b(calibrationTool,level1);
end

disp('Integration successful')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

end
