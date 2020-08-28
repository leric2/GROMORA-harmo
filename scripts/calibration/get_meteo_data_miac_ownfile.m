function meteo = get_meteo_data_miac_ownfile(calibrationTool)
%==========================================================================
% NAME          | get_meteo_data_miac_ownfile
% TYPE          | function
% AUTHOR(S)     | Adapted from 'read_meteo_data_ownfile.m', by F. Schranz
% CREATION      | 2020-06-12
%               |
% ABSTRACT      | Function to read meteo data from the MIAWARA-C weather
%               | station when it is saved in a separate file
%               | MIAWARA-C_meteo*.txt
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: calibrationTool
%               |
%               | OUTPUTS: meteo
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
%%

file=calibrationTool.meteoFile;

% initialize return value
meteo.file = file;

files = dir([calibrationTool.meteoFile '.txt']);

meteo.time = [];
meteo.dateTime = [];
meteo.temperature=[];
meteo.air_temperature = [];
meteo.rel_humidity= [];
meteo.pressure=[];
meteo.wind_speed= [];
meteo.wind_dirn = [];
meteo.visibility= [];
meteo.global_rad =  [];
meteo.precipitation=  [];
	
if ~isempty(files) 
	
    [y, ~] = readtext([file '.txt'], calibrationTool.delimiter_logfile, '', '"');
    if ischar(cell2mat(y(end,1))), y=y(1:end-1,:); end %% delete last row (power break)

    header = y(1,:);
    u = length(header);

    for k = 1:15 %column 16: string, 	
        eval( sprintf('%s_1= cell2mat(y(2:end,k));',  char(header(k))));
        eval( sprintf('log1.%s  = cell2mat(y(2:end,k));',  char(header(k))) );
    end

    date_1 = '';
    for k = 1:length(Month_1)
        date_1 = [date_1;datestr(datenum(sprintf('%i/%i/%i\n', Month_1(k),Day_1(k),Year_1(k)))+(Hour_1(k)+ Minute_1(k)/60+Second_1(k)/3600)/24,31)];
    end

    log1.time = date_1;
        
	meteo.time_str      = log1.time;
    meteo.time            = datenum(meteo.time_str);
    meteo.dateTime      = meteo.time;
	meteo.temperature   = log1.Meteo.TAmb_C+273.15;
    meteo.air_temperature = meteo.temperature;
	meteo.rel_humidity  = log1.Meteo.RH;
	meteo.pressure      = log1.Meteo.Pressure_hPA;
    meteo.air_pressure  = meteo.pressure;
	meteo.wind_speed    = log1.Wind.Current;
	meteo.wind_dirn     = log1.Wind.DirN;
	meteo.visibility    = 0*ones(1,length(meteo.time));
	meteo.global_rad    = 0*ones(1,length(meteo.time));
	meteo.precipitation = log1.Meteo.Precipitation_mm;
	     
end
