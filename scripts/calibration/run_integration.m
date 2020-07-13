function calibrationTool = run_integration(calibrationTool)
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
% Defining level1a filename to read (to be adapted for other users)
filename = [calibrationTool.level1Folder calibrationTool.instrumentName '_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];
calibrationTool.filenameLevel1a = filename;

% Check that the calibration file (level1a) exist for this day
if isfield(calibrationTool,'filenameLevel1a')
    if ~exist(calibrationTool.filenameLevel1a,'file')
        error('No calibration data found for this day')
    end
else
    error('No calibration data found for this day')
end

% Start integration
disp(['Starting the integration process for ' calibrationTool.instrumentName ': ' calibrationTool.dateStr])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1a to level 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We create level1b structure which will contain both the calibrated
% spectra read from the level1a data (level1b.calibratedSpectra and the 
% integrated spectra that will be added later (level1b.integration) 
level1b = struct();

[level1b.calibratedSpectra, meteoData, calibrationTool] = calibrationTool.read_level1a(calibrationTool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now because no Payerne dataset
% calibrationTool.meteoFolder = '/mnt/instrumentdata/meteo/exwi/meteo/';
% calibrationTool.get_meteo_data  =  @(calibrationTool,correctedSpectra) get_meteo_data_unibe(calibrationTool,correctedSpectra);

if isempty(meteoData)
    %TODO
    disp('no meteo data for this day...')
else
    level1b.calibratedSpectra = calibrationTool.add_meteo_data(calibrationTool, meteoData, level1b.calibratedSpectra);
end

% checking the quality of the channels and flagging the potential bad ones
% (we do not remove any)
level1b.calibratedSpectra = calibrationTool.checking_channel_quality(level1b.calibratedSpectra,calibrationTool,1);

% Compute tropospheric transmittance and correction for every calibrated
% spectra.
level1b.calibratedSpectra = calibrationTool.tropospheric_correction(level1b.calibratedSpectra, calibrationTool,10.4);

% Integrating the "good spectra" based on tropospheric transmittance and
% calibration flags. --> To improve. Maybe introduce weighted mean of
% spectra based on tropospheric transmittance ?
level1b.integration = calibrationTool.integrate_calibrated_spectra(calibrationTool,level1b.calibratedSpectra);

%% Correction and checks
% Now on the integrated spectra; checking the quality of the channels and 
% flagging the potential bad ones (we do not remove any).
level1b.integration = calibrationTool.checking_channel_quality(level1b.integration,calibrationTool,2);

% Performing window correction
level1b.integration = calibrationTool.window_correction(calibrationTool,level1b.integration);

% Compute tropospheric transmittance and correction for every integrated
% spectra.
level1b.integration = calibrationTool.tropospheric_correction(level1b.integration, calibrationTool,10.4);

% Check integrated spectra and define the flags for level1b
level1b.integration = calibrationTool.check_integrated(calibrationTool, level1b.integration);

% sideband correction ?
% TODO

% Plotting and saving calibrated and corrected spectra
if calibrationTool.integratedSpectraPlot
    calibrationTool.plot_integrated_spectra(calibrationTool,level1b.integration,50,260)
end

%%
% Saving integrated spectra (level1b) into NetCDF-4 file
disp('Saving Level 1b...')
if calibrationTool.numberOfSpectrometer > 1
    %TODO save the 4 integrated spectra in 1 level1b
else
    calibrationTool  =  calibrationTool.save_level1b(calibrationTool,level1b);
end

disp('Integration successful')
calibrationTool.successfulIntegration = true;

end
