function y = iap_read_2019(file, hk_only)

% y = iap_read(file, hk_only)
%
% Reads houskeeping and binary data of the universal IAP data format
%
% INPUT
% file:     Filename without extension
%
% OUTPUT
% y:       Housekeeping structure with field names for each parameter
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
% 2017-??-??    A. Murk (adapted for MOPI5)
% 2019-01-03    A. Murk (adapted for new MOPI5 LV version)

if nargin<2 hk_only=0; end

% initialize return value
y.file = file;
y.comment = [];

if findstr(file, '.txt')
    file=file(1:end-4);
end

% ========= read logfile ===================

% x = readtext([file '.txt'], ';');          % function found under www.mathworks.nl/matlabcentral/fileexchange/10946-readtext
x   = importdata([file '.txt'], ';', 1); % mucht faster than readtext if HK data is only numbers

[M,N] = size(x.data);  % M = number of spectra, N = number of columns

for i = [1:N]
    name = x.colheaders{i};
    name(name=='.' | name==' ') = '_'; % remove blanks and dots from field names
    value = x.data(:,i); 
    if value==value(1)
        value=value(1);
    end  
    y = setfield( y, name, value );
end

for i=0:7
    name  = sprintf('NrOfChannels_%d',i); 
    if isfield(y, name)
       value = getfield(y,name); 
       if length(value)==1;
          y.NrOfChannels(i+1)=value;
       else
           y.NrOfChannels(i+1)=NaN;
           warning('Missing Spectra!')
       end
    end
end


% disp(y.NrOfChannels);


% % susanna
% date = '';
% for k = 1:M
%     date = [date;datestr(datenum(sprintf('%i/%i/%i\n', log.Month(k),log.Day(k),log.Year(k)))+(log.Hour(k)+ log.Minute(k)/60+log.Second(k)/3600)/24,31)];
% end
% log.date = date;
% log.time = datenum(date);
% save([name2 '_log' ] ,'log');

% calculate time in [h]
if isfield(y, {'Hour' 'Minute' 'Second'})
    y.t = y.Hour + y.Minute/60 + y.Second/3600;
end

if hk_only
    y.data = NaN;
else
    
    D = dir([file '.bin']);
    
    if isempty(D)
        y.data=NaN;
        return
    else
         channels= sum(y.NrOfChannels);   
%         D.bytes/4/channels; % 4 bytes for each floating point value
    end
    
    % read complete binary data in one array (huge!)
    fid = fopen( [file '.bin'], 'r', 'b');
    y.data = fread(fid, [channels,M], 'float32');
    fclose(fid);    
end

% settings for MOPI tests 
% 2017
% y.T = [y.HousingParameters_Temp_5, y.HousingParameters_Temp_6, y.HousingParameters_Temp_7];
% 2019
y.T = [y.HousingParameters_Temp_0, y.HousingParameters_Temp_7]; 



