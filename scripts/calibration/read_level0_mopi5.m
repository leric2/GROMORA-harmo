function [log, data] = read_level0_mopi5(calibrationTool, rawFileReading)
%==========================================================================
% NAME          | read_level0_generic.m
% TYPE          | function
% AUTHOR(S)     | Axel Murk (adapted to GROSOM by E.S.)
% CREATION      | 2011-03-04 (adapted 01.2020)
%               |
% ABSTRACT      | Reads houskeeping and binary data of the universal IAP 
%               | data format
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: Filename without extension
%               |
%               | OUTPUTS:  1. Housekeeping structure with field names for each parameter
%               |           2. [length(file)*#channels] line vector of binary file
%               |             
% CALLED by     | run_retrieval(calibrationTool);
%               |
% CALLS         | readtext([file '.txt'], calibrationTool.delimiter_logfile, '', '"');
%               |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FORMAT:
%
% file.txt:
% -  Arbitrary number of lines with comments starting with '%'
%    (e.g. instrument, software version, commnets of the observer...)
% -  One line with the name of M housekeeping parameters separated by ' '
% -  N lines of the houskeeping values with M entries each
% -  Adapted for files separated with given delimiter e.g. ';' (F. Schranz 2020-06-08)
%
% file.bin
% 32bit floatinig point data 
%==========================================================================
if nargin<2 rawFileReading=0; end

file=calibrationTool.file;
% initialize return value
log.file = file;
log.comment = [];

% ========= read logfile ===================
fid = fopen( [file calibrationTool.logFileDataExtension], 'r');

while 1
    s = fgetl(fid);
    if findstr('%%', s)
        log.comment=[log.comment '\n' s];
    else
        break
    end
end

if s(1)=='%';  s(1)=[]; end

if isfield(calibrationTool,'delimiter_logfile')
    %header = textscan(s, '%s','delimiter', calibrationTool.delimiter_logfile);
    %header = header{1}; % cell array with all header parameters
    %N = length(header); % number of header parameters

    %[y, result] = readtext([file '.txt'], calibrationTool.delimiter_logfile, '', '"');
    y = importdata([file calibrationTool.logFileDataExtension], calibrationTool.delimiter_logfile, 1);
    header = y.colheaders;
    N = length(header);
    x = y.data;
    %if ischar(cell2mat(y(end,1))), y=y(1:end-1,:); end

    %x = cell2mat(y(2:end,:))';
    
else
    error('Please provide the delimiter symbol for the log file')
end

log.x = x; 
log.header = header; 

M = size(x,1);     % number of data entries
%fclose(fid);

for n = 1:N
    name = header{n}; 
    name(name=='.')='_'; 
    name(name==' ')='_'; 
    if length(name)<2 continue; end; 
    %disp(name)
    log = setfield(log, name, x(:,n));
end

% calculate time in [h]
if isfield(log, {'Hour' 'Minute' 'Second'})
    log.t = log.Hour + log.Minute/60 + log.Second/3600;
end

% calculate time in [s] ?? also hours?
if isfield(log, {'Hour' 'Min' 'Sec'})
    log.t = log.Hour + log.Min/60 + log.Sec/3600;
end

if rawFileReading
    K=[1 2 5 7]; % Generic MOPI software includes 8 optional FFTS. Only these 4 are used
    ffts_model = calibrationTool.ffts_model;
    
    if ffts_model==0 return; end
    
for i = ffts_model
		
    % vector with file pointers to the beginning of each spectrum
    bos = getfield( log, sprintf('BinFilePointer_%d', K(i)) );
    N   = length(bos); % number of cycles
    
    CHANNELS = getfield( log, sprintf('NrOfChannels_%d', K(i)) );  % Number of channels
    if length(CHANNELS)==1
        channels=CHANNELS;
        missing_channels=[];
    else
        channels=max(CHANNELS);
        missing_channels=find(CHANNELS<channels);
    end
    
    
    data = single(zeros(N,channels));
    fid = fopen( [calibrationTool.file '.bin'], 'r', 'b');
    
    h = waitbar(0,sprintf('Reading %s in %s.bin', calibrationTool.spectrometer, file));
    for n=1:N
        fseek(fid, bos(n),-1);
        data(n,:) = fread(fid, channels, 'float32');
        waitbar(n/N,h)
    end
    close (h);
    
    % mark cycles with missing channels as NaN
    data(missing_channels, :)=NaN;
    fclose(fid);
    % 	% generate Frequency axis
    % 	if i<3
    % 		log.f = FS(i)/2 * [-1:2/channels:1-2/channels]; % I/Q processing in USRP
    % 	else
    % 		log.f = FS(i)/2 * [0:1/channels:1-1/channels]; % real data in AC240 and U5303
    % 	end
    
end
else
    data = NaN;
end

% %% read missing measurements which were saved in the previous file
% 
% file2 = [calibrationTool.extraFileFolder calibrationTool.instrumentName '_' datestr(datenum(calibrationTool.dateStr)-1,'YYYYmmdd') '_extra'];
%     
% if exist(file2,'file')
%     
%     
%     log2.file = file2;
%     log2.comment = [];
% 
% 
%     % ========= read logfile ===================
%     fid = fopen( [file2 '.txt'], 'r');
% 
%     while 1
%         s = fgetl(fid);
%         if findstr('%%', s)
%             log2.comment=[log2.comment '\n' s];
%         else
%             break
%         end
%     end
% 
%     if s(1)=='%';  s(1)=[]; end
% 
%     if isfield(calibrationTool,'delimiter_logfile')
%         header = textscan(s, '%s','delimiter', calibrationTool.delimiter_logfile);
%         header = header{1}; % cell array with all header parameters
%         N = length(header); % number of header parameters
% 
%         [y, result] = readtext([file2 '.txt'], calibrationTool.delimiter_logfile, '', '"');
%         if ischar(cell2mat(y(end,1))), y=y(1:end-1,:); end
% 
%         x = cell2mat(y(2:end,:))';
% 
%     else
%         header = textscan(s, '%s'); 
% 
%         header = header{1}; % cell array with all header parameters
%         N = length(header); % number of header parameters
%         % x = fscanf(fid, '%f;', [N, inf]);  % data array
% 
%         x = fscanf(fid, '%f ', [N, inf]);  % data array
% 
%     end
% 
% 
%     M = size(x,2);     % number of data entries
%     fclose(fid);
% 
%     for n = 1:N
%         name = header{n}; 
%         name(name=='.')='_'; 
%         name(name==' ')='_'; 
%         if length(name)<2 
%             continue; 
%         end
%         %disp(name)
%         log2 = setfield(log2, name, x(n,:));
%     end
% 
%     % calculate time in [h]
%     if isfield(log2, {'Hour' 'Minute' 'Second'})
%         log2.t = log2.Hour + log2.Minute/60 + log2.Second/3600;
%     end
% 
%     % calculate time in [s] ?? also hours?
%     if isfield(log2, {'Hour' 'Min' 'Sec'})
%         log2.t = log2.Hour + log2.Min/60 + log2.Sec/3600;
%     end
% 
% 
%     D = dir([file2 '.bin']);
% 
%     if isempty(D)
%         rawSpectra_pre=NaN;
%     else
%         theoreticalNumberDataEntries=D.bytes/calibrationTool.numberOfChannels/4;    % 4 bytes for each floating point value
%     end
% 
%     % read complete binary data in one vector
%     if nargout>1
%         fid = fopen( [file2 '.bin'], 'r', calibrationTool.binaryType);
%         rawSpectra_pre = fread(fid ,calibrationTool.numberOfChannels*theoreticalNumberDataEntries,'float32=>float32');
%         fclose(fid);
%     end
% 
%     % we want a line vector for the following
%     rawSpectra_pre=rawSpectra_pre';
% 
%     log2.x = x; 
%     log2.header = header; 
% 
% 
%     % concatenate
% 
%     rawSpectra_out = [rawSpectra_pre rawSpectra];
%     logFile_out    = cell2struct(cellfun(@vertcat,struct2cell(log2),struct2cell(log),'uni',0),fieldnames(log),1);
% 
% else
%     rawSpectra_out = rawSpectra;
%     logFile_out    = log;
% end






