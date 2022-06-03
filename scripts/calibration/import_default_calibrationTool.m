function calibrationTool = import_default_calibrationTool(dateStr)
%==========================================================================
% NAME          | import_default_calibrationTool(instrumentName)
% TYPE          | Function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Creation of the calibrationTool structure for GROSOM
%               | project. At this point, we just initiate the structure
%               | with some basic parameters (dates and physical constant).
%               | The complete definition differs for each instrument and 
%               | is done in "import_instrumentName_calibrationTool.m". 
%               | 
%               |
% ARGUMENTS     | INPUTS:  1. dateStr: in the form YYYY_MM_DD
%               |
%               | OUTPUTS: 1. calibrationTool: the default toolbox for
%               |               launching a calibration or integration for
%               |               for this day.
%               |
%==========================================================================
calibrationTool = struct();

% Date and time parameters
calibrationTool.dateStr=dateStr;

calibrationTool.Year = str2double(dateStr(1:4));
calibrationTool.Month = str2double(dateStr(6:7));
calibrationTool.Day = str2double(dateStr(9:10));

% Conversion of datestr to datenum:
calibrationTool.dateTime = datetime(dateStr,'InputFormat','yyyy_MM_dd');
calibrationTool.timeNumber=datenum(dateStr,'yyyy_mm_dd');

% Name of the instrument
calibrationTool.dateStr=dateStr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constant
calibrationTool.lightSpeed=299792458;   % [m/s] 
calibrationTool.h = 6.62606876e-34;       % [Js]
calibrationTool.kb=1.38065e-23;         % [J/K]

% Conversion between degree C and Kelvin.
calibrationTool.zeroDegInKelvin = 273.15;

%calibrationTool.deltaTCorr = 10.4;

% MW background radiation at TOA
calibrationTool.backgroundMWTb = 2.736;
end
