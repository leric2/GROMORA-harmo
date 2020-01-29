function [x,data]=mopi5_read(retrievalTool)
%function x=mopi5_read(file, ffts_model)
tic

file=retrievalTool.file;
ffts_model=retrievalTool.ffts_model;

S  = {'USRP-A', 'USRP-B','U5303', 'AC240'}; % FFTS model 1 to 4
FS = [200 20  3200 2000]; % sampling rates in MHz 

%x=iap_read_2019(file,1);

% initialize return value
clear x
x.file = file;
x.comment = [];

% ========= read logfile ===================
fid = fopen( [file '.txt'], 'r');

while 1
    s = fgetl(fid);
    if findstr('%%', s)
        x.comment=[x.comment '\n' s];
    else
        break
    end
end

if s(1)=='%';  s(1)=[]; end

header = textscan(s, '%s','delimiter', ';');
%header = textscan(s, '%s'); 
header = header{1}; % cell array with all header parameters
N = length(header); % number of header parameters
% x = fscanf(fid, '%f;', [N, inf]);  % data array
log = fscanf(fid, '%f; ', [N, inf]);  % data array
M = size(x,2);     % number of data entries
fclose(fid);

for n = 1:N
    name = header{n}; 
    name(name=='.')='_'; 
    name(name==' ')='_'; 
    %if length(name)<2 continue; end; 
    %disp(name)
    x = setfield(x, name, log(n,:));
end

% calculate time in [h]
if isfield(x, {'Hour' 'Minute' 'Second'})
    x.t = x.Hour + x.Minute/60 + x.Second/3600;
end

% calculate time in [s]
if isfield(x, {'Hour' 'Min' 'Sec'})
    x.t = x.Hour + x.Min/60 + x.Sec/3600;
end

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
	
	
	data = zeros(N,channels);
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
toc








