function rawSpectra = reformat_spectra_generic(rawSpectra,logFile,calibrationTool)
    % From vector to matrix
    % Number of spectra for this day
    %n=size(log.t,1);
    n=length(logFile.t);
    % reshaping the initial raw vector
    try
        rawSpectra=(reshape(rawSpectra,[calibrationTool.numberOfChannels,n]))';
    catch ME
        error(ME.identifier,'Problem when reformatting the spectra');
        disp(ME.message)
    end
end

