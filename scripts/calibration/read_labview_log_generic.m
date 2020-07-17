function labviewLog = read_labview_log_generic(instrumentName)
%==========================================================================
% NAME          | read_labview_log_generic
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
%               |
% ABSTRACT      | Function to read
%               | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:
%               |
%               | OUTPUTS:
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================

labviewLogFolder = '/home/esauvageat/Documents/GROSOM/Analysis/InputsCalibration/';
labviewLogFile=[labviewLogFolder instrumentName '_Labview_Log.txt'];
    
% Transforming it into matlab structure
fileID = fopen(labviewLogFile,'r');
f = textscan(fileID,'%s','Delimiter','\r', 'EndOfLine','\n');

%log = cell2struct(f,'te') 

for i = 1:length(f{1, 1})
    line = f{1, 1}(i);
    if ~isempty(line{1,1})
        try
            dateN =  datenum(line{1,1}(1:19),'dd.mm.yyyy HH:MM:SS');
            labviewLog(i).dateNum = dateN-datenum(1970,1,1);
            labviewLog(i).dateTime = datetime(dateN ,'ConvertFrom','datenum');
            labviewLog(i).message = line{1,1}(20:end);
        catch ME
            warning(['Problem with log reading line ' num2str(i)])
        end
    end
end