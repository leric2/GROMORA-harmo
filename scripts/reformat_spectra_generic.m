function rawSpectra = reformat_spectra_generic(rawSpectra,log,retrievalTool)
% From vector to matrix
% Number of spectra for this day
n=size(log.x,2);
% reshaping the initial raw vector
rawSpectra=(reshape(rawSpectra,[retrievalTool.numberOfChannels,n]))';
end

