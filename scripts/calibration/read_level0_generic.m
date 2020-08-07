function [log,rawSpectra] = read_level0_generic(calibrationTool)
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
% CALLED by     | run_retrieval(retrievalTool);
%               |
% CALLS         | readtext([file '.txt'], retrievalTool.delimiter_logfile, '', '"');
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
    header = textscan(s, '%s','delimiter', calibrationTool.delimiter_logfile);
    header = header{1}; % cell array with all header parameters
    N = length(header); % number of header parameters

    [y, result] = readtext([file calibrationTool.logFileDataExtension], calibrationTool.delimiter_logfile, '', '"');
    
    if ischar(cell2mat(y(end,1))), y=y(1:end-1,:); end
    
    if calibrationTool.positionIndAsName
        indexPos = find(strcmp(header,'Position'));
        newPositionColumn = nan*ones(length(y(2:end,27)),1);
        newPositionColumn(strcmp(y(2:end,27),calibrationTool.nameColdIndice)) = calibrationTool.indiceCold;
        newPositionColumn(strcmp(y(2:end,27),calibrationTool.nameHotIndice)) = calibrationTool.indiceHot;
        newPositionColumn(strcmp(y(2:end,27),calibrationTool.nameAntennaIndice)) = calibrationTool.indiceAntenna;
        newPositionColumn(strcmp(y(2:end,27),calibrationTool.otherName)) = -1;
        y(2:end,27) = num2cell(newPositionColumn);
    end
    x = cell2mat(y(2:end,:))';
else
    header = textscan(s, '%s'); 
    
    header = header{1}; % cell array with all header parameters
    N = length(header); % number of header parameters
    % x = fscanf(fid, '%f;', [N, inf]);  % data array

    x = fscanf(fid, '%f ', [N, inf]);  % data array
    
end


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

% calculate time in [s] ?? also hours?
if isfield(log, {'Hour' 'Min' 'Sec'})
    log.t = log.Hour + log.Min/60 + log.Sec/3600;
end

D = dir([file calibrationTool.binaryDataExtension]);

if isempty(D)
    rawSpectra=NaN;
else
    theoreticalNumberDataEntries=D.bytes/calibrationTool.numberOfChannels/4;    % 4 bytes for each floating point value
end

% read complete binary data in one vector
if nargout>1
    fid = fopen( [file calibrationTool.binaryDataExtension], 'r', calibrationTool.binaryType);
    rawSpectra = fread(fid,calibrationTool.numberOfChannels*theoreticalNumberDataEntries,'float32=>float32');
    fclose(fid);
end

% we want a line vector for the following
if (size(rawSpectra,2)==1)
    rawSpectra=rawSpectra';
end

log.x = x; 
log.header = header; 
end

