function rawSpectra = flip_spectra_gromos(rawSpectra)
% frequency order in the IQ spectra is [0:-500MHz and 0:500]
% flip the negative frequencies (the first 16384 channels) to obtain -500:500MHz

midChannel=16384;
rawSpectra = [fliplr(rawSpectra(:,1:midChannel)) rawSpectra(:,midChannel+1:end)];
end

