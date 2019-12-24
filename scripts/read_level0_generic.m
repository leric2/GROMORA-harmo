function [log,rawSpectra] = read_level0_generic(file,retrievalTool)
% [log,rawSpectra,readingLevel0Error] = read_level0_generic(file)
%
% Reads houskeeping and binary data of the universal IAP data format
%
% INPUT
% file:     Filename without extension
%
% OUTPUT
% hk:       Housekeeping structure with field names for each parameter
% data:     [N,channles] array of binary file
%
%
% DATA FORMAT
%
% file.txt:
% -  Arbitrary number of lines with comments starting with '%'
%    (e.g. instrument, software version, commnets of the observer...)
% -  One line with the name of M housekeeping parameters separated by ' '
% -  N lines of the houskeeping values with M entries each
%
% file.bin
% 32bit floatinig point data for an array of [Nxchannles] values
%
% The number of channles is determined from file size

% Revisions
% 2011-03-04    A. Murk (first issue)


% initialize return value
log.file = file;
log.comment = [];

% ========= read logfile ===================
fid = fopen( [file '.txt'], 'r');

while 1
    s = fgetl(fid);
    if findstr('%%', s)
        log.comment=[log.comment '\n' s];
    else
        break
    end
end

if s(1)=='%';  s(1)=[]; end

% header = textscan(s, '%s','delimiter', ';');
header = textscan(s, '%s'); 
header = header{1}; % cell array with all header parameters
N = length(header); % number of header parameters
% x = fscanf(fid, '%f;', [N, inf]);  % data array
x = fscanf(fid, '%f ', [N, inf]);  % data array
M = size(x,2);     % number of data entries
fclose(fid);

for n = 1:N
    name = header{n}; 
    name(name=='.')='_'; 
    name(name==' ')='_'; 
    if length(name)<2 continue; end; 
    %disp(name)
    log = setfield(log, name, x(n,:));
end

% calculate time in [h]
if isfield(log, {'Hour' 'Minute' 'Second'})
    log.t = log.Hour + log.Minute/60 + log.Second/3600;
end

% calculate time in [s]
if isfield(log, {'Hour' 'Min' 'Sec'})
    log.t = log.Hour + log.Min/60 + log.Sec/3600;
end


D = dir([file '.bin']);

if isempty(D)
    rawSpectra=NaN;
else
    theoreticalNumberDataEntries=D.bytes/retrievalTool.numberOfChannels/4;    % 4 bytes for each floating point value
end

% read complete binary data in one vector
if nargout>1
    fid = fopen( [file '.bin'], 'r', retrievalTool.binaryType);
    rawSpectra = fread(fid ,retrievalTool.numberOfChannels*theoreticalNumberDataEntries,'float32=>float32');
    fclose(fid);
end

% we want a line vector for the following
rawSpectra=rawSpectra';

log.x = x; 
log.header = header; 
end

