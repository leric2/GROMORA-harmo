
function write_tau_and_Teff(date,instrumentName)

dateStr=datestr(date,'yyyy_mm_dd');
    
% Import default tools for running a retrieval for a given instrument
calibrationTool=import_default_calibrationTool(instrumentName,dateStr);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading and formatting the raw spectra for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading level0 data...')

% Reading raw data
[logFile,rawSpectra] = calibrationTool.read_level0(calibrationTool,1);



% The raw log file from each instrument is different and we should try to
% harmonize it as much as possible (different function for each
% instrument and might need date information later ?).
logFile = calibrationTool.harmonize_log(calibrationTool,logFile);

% Quality check of the raw data:
if calibrationTool.checkLevel0
    warningLevel0 = calibrationTool.check_level0(logFile,rawSpectra,calibrationTool);
else
    warningLevel0 = '';
end

%% TO CHECK IF RIGHT
% Reformat the raw spectra from vector to matrix
if size(rawSpectra,1) == 1
    rawSpectra = calibrationTool.reformat_spectra(rawSpectra,logFile,calibrationTool);
end

% when needed, flip it !
if calibrationTool.flipped_spectra
    rawSpectra = calibrationTool.flip_spectra(rawSpectra);
end

% TODO
% remove channels which are known to be bad (2pol)

% TODO
% check_level0


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_meteo_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logFile.meteo = calibrationTool.read_meteo_data(calibrationTool);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tipping Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logFile.TC = calibrationTool.run_tipping_curve(rawSpectra, logFile, calibrationTool);
