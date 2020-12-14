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
%           | calibrate. 
%           | 
%           |
% ARGUMENTS | INPUTS: - calibrationTool: structure containing all
%           | information about the integration we want to perform.
%           | Documentation about this structure can be found in external
%           | document.
%           |
%           |
%           | OUTPUTS: - Level1b
%           |           
%           |
%           |
% CALLS     | Some depends on instruments, all are stored in calibrationTool:
%           | %%%%%%%%%%%%%%%%%%%%% Level1a -> Level1b
%           | read_level1a
%           | get_meteo_data
%           | checking_channel_quality
%           | tropospheric_correction_generic
%           | window_correction
%           | plot_integrated_spectra
%           | save_level1b
%           |
%==========================================================================
% Check that the calibration file (level1a) exist for this day
if ~exist(calibrationTool.filenameLevel1a,'file')
    error('No calibration data found for this day')
 end

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

if isempty(meteoData)
    %TODO
    disp('no meteo data for this day...')
else
    level1.calibratedSpectra = calibrationTool.add_meteo_data(calibrationTool, meteoData, level1.calibratedSpectra);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction and integration of the calibrated spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking the quality of the channels and flagging the potential bad ones
% (we do not remove any)
filterType = calibrationTool.filterTypeChannelQualityCal;
level1.calibratedSpectra = calibrationTool.checking_channel_quality(level1.calibratedSpectra,calibrationTool,filterType);

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
% Now on the integrated spectra; checking the quality of the channels and 
% flagging the potential bad ones (we do not remove any).
filterType = calibrationTool.filterTypeChannelQualityInt;
level1.integratedSpectra = calibrationTool.checking_channel_quality(level1.integratedSpectra,calibrationTool,filterType);

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
if calibrationTool.integratedSpectraPlot && nansum([level1.integratedSpectra.outlierCalib])<length(level1.integratedSpectra)
    calibrationTool.plot_integrated_spectra(calibrationTool,level1.integratedSpectra)
elseif calibrationTool.integratedSpectraPlot && length(level1.integratedSpectra)==1
    calibrationTool.plot_integrated_spectra(calibrationTool,level1.integratedSpectra)
end

% Saving integrated spectra (level1b) into NetCDF-4 file
disp('Saving Level 1b...')
calibrationTool  =  calibrationTool.save_level1b(calibrationTool,level1);

disp('Integration successful')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

end
