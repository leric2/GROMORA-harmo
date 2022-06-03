function [logFile, rawSpectra] = read_level0_wirac(calibrationTool, readRawFile)

% y = iap_read(file, hk_only, channels, index)
%
% Reads houskeeping (*.txt) and binary (*.bin) data of the IAP radiometers
%
% INPUT
% file:     Filename (without  *.txt or *.bin)
% hk_only:  Flag to read only houskeeping data (default 0)
% channels: Number of spectrometer channels in the *.bin file
%           default 0, in this case number of channels is determined from file size
% index:    Read the entire file (index=0),
%           only the first 1:index measurements (if inedx has length 1), or
%           the requested range of measurements if length(index)>1
%
% OUTPUT
% y        Structure with field names for each houskeeping parameter
% y.x      [N,channles] array with data from the binary file (NaN if hk_only=0)
%
%
% DATA FORMAT
%
% file.txt:
% -  Arbitrary number of lines with comments starting with
%    (e.g. instrument settings, software version, date+time, ...)
% -  One line with the name of M housekeeping parameters
% -  N lines of the houskeeping values with M entries each
%
% file.bin
% 32bit floatinig point data for an array of [Nxchannels] values
%
% The number of channels is determined from file size

% Revisions
% 2011-03    A.M. (first issue)
% 2021-03    A.M. Updated with readtable to be more flexible with field separators
% 2021-07    A.M. Updated to check whether *.bin file can be read in one or line by line in case the data is corrupted
% 
% if nargin<2    hk_only =0; end % 1: read only housekeeping data file
% if nargin<3    channels=0; end % 0: determin number of channels from file size
% if nargin<4    index =0; end % 0: read all measurements

file = calibrationTool.file;
channels = calibrationTool.numberOfChannels;
% make sure that file name is given without exetension
if findstr(file, '.txt');
    file=file(1:end-4);
end

% initialize return value
y.file = file;


% read housekeeping data and determine nuber of measurements M
x=readtable(file);

[M,N]=size(x);

% determine number of spectrometer channels from file size
D = dir([file '.bin']);
if isempty(D)
    hk_only=1;
    clear D; 
    D.bytes=NaN;
end
LineByLine=0;
if channels > 0
    % check if *.bin file size matches to given number of channels and expected number of measurements
    if D.bytes ~= channels*4*M
        warning( 'Size of binary file %d does not match channels*measurements*4=%d\nDAta wil be read line by line', D.bytes, channels*4*M)
        LineByLine=1;
    end
else
    % estimate number of channels from fie size and M
    channels=D.bytes/4 /M; % 4 bytes for each floating point value
end



% use only the first entries up to index
% if index
%     if length(index)==1
%         index=1:index;
%     end    
%     x=x(index,:);
%     M=length(index);
% end

% fprintf('channels= %d    Spectra= %d\n', channels, M);



% convert table to structures
for n = 1:N
    name = x.Properties.VariableNames{n};
    data = x{:,n};
    if iscell(data)
         % data=data{1};
    else
        if max(diff(data))==0
            data=data(1);
        end
    end
    y = setfield(y, name, data);
end


% calculate time in [h]
if isfield(y, {'Hour' 'Minute' 'Second'})
    y.t = y.Hour + y.Minute/60 + y.Second/3600;
end

% calculate time in [h]
if isfield(y, {'Hour' 'Min' 'Sec'})
    t = y.Hour + y.Min/60 + y.Sec/3600;
    k=find(diff(t)<0);
    for i=1:length(k)
        t( k(i)+1:end ) = t( k(i)+1:end )+24;
    end
    y.t = t;
end


if isfield(y, 'DateTime')
    DN     = datenum(x.DateTime, 'yyyy-mm-ddTHH:MM:SS.FFF');    
    y.DN   = DN; 
    y.date = datestr(DN(1), 'yyyy-mm-dd');
    t0     = DN(1); 
    t0     = hour(t0) + minute(t0)/60 + second(t0)/3600;
    y.t    = (DN-DN(1))*24 + t0;     
end


if ~readRawFile
    y.x=NaN;
    return
end

if isfield(y, 'size_before')      pos=y.size_before;      end
if isfield(y, 'BinFilePointer_0') pos=y.BinFilePointer_0; end



y.x = NaN(channels,M, 'single');


fid = fopen( [file '.bin'], 'r', 'b');


% set file pointer to the first entry which shall be read, only needed if a subset of the file shall be read
% if length(index)>1
%     fseek(fid, pos(1), 'bof');
% end


if 0
    h = waitbar(0,sprintf('Reading %s in %s.bin', M, file));
	for i=1:M
        waitbar(i/M,h)        
        fseek(fid, pos(i), 'bof');
        y.x(:,i) = fread(fid, channels, 'float32=>single');
    end
    close(h);   
else
    y.x = fread(fid, [channels,M], 'float32=>single');  
end
fclose(fid);
logFile = y;
rawSpectra = y.x';
end



