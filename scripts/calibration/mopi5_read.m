function [x,data]=mopi5_read(calibrationTool)
%function x=mopi5_read(file, ffts_model)
tic

file=calibrationTool.file;
ffts_model=calibrationTool.ffts_model;

S  = {'USRP-A', 'USRP-B','U5303', 'AC240'}; % FFTS model 1 to 4
FS = [200 20  3200 2000]; % sampling rates in MHz 

% reading the log file only:
x=iap_read_2019(file,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=[1 2 5 7]; % Generic MOPI software includes 8 optional FFTS. Only these 4 are used

if ffts_model==0 return; end

for i = ffts_model
		
	% vector with file pointers to the beginning of each spectrum
	bos = getfield( x, sprintf('BinFilePointer_%d', K(i)) ); 
	N   = length(bos); % number of cycles
	
	CHANNELS = getfield( x, sprintf('NrOfChannels_%d', K(i)) );  % Number of channels
	if length(CHANNELS)==1
		channels=CHANNELS;
		missing_channels=[]; 
	else
		channels=max(CHANNELS);
		missing_channels=find(CHANNELS<channels);
	end
	
	
	data = single(zeros(N,channels));
	fid = fopen( [file '.bin'], 'r', 'b');
	
	h = waitbar(0,sprintf('Reading %s in %s.bin', S{i}, file));
	for n=1:N
 		fseek(fid, bos(n),-1);
		data(n,:) = fread(fid, channels, 'float32');
		waitbar(n/N,h)
	end
	close (h);
	
	% mark cycles with missing channels as NaN
	data(missing_channels, :)=NaN; 
	
	% generate Frequency axis 
	if i<3
		x.f = FS(i)/2 * [-1:2/channels:1-2/channels]; % I/Q processing in USRP
	else
		x.f = FS(i)/2 * [0:1/channels:1-1/channels]; % real data in AC240 and U5303
	end
	
end
fclose(fid);
%data=tall(data);
toc








