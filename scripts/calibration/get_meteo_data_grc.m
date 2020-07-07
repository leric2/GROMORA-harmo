function meteo = get_meteo_data_grc(calibrationTool)
%==========================================================================
% NAME          | get_meteo_data_grc
% TYPE          | function
% AUTHOR(S)     | Adapted from 'read_meteo_data_gromosc', by F. Schranz
% CREATION      | 2020-06-12
%               |
% ABSTRACT      | Function to read meteo data from the GROMOS-C weather
%               | station
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
meteo.temperature=[];
meteo.rel_humidity= [];
meteo.pressure=[];
meteo.wind_speed= [];
meteo.wind_dirn = [];
meteo.visibility= [];
meteo.global_rad =  [];
meteo.precipitation=  [];
	
if ~isempty(files) 
	
    [y, ~] = readtext([files.folder '/' files.name], calibrationTool.delimiter_logfile, '', '"');
    if ischar(cell2mat(y(end,1))), y=y(1:end-1,:); end %% delete last row (power break)

    header = y(1,:);
    u = length(header);

    for k = 1:u %column 16: string, 	
        eval( sprintf('%s_1= cell2mat(y(2:end,k));',  char(header(k))));
        eval( sprintf('log1.%s  = cell2mat(y(2:end,k));',  char(header(k))) );
    end

    date_1 = '';
    for k = 1:length(Month_1)
        date_1 = [date_1;datestr(datenum(sprintf('%i/%i/%i\n', Month_1(k),Day_1(k),Year_1(k)))+(Hour_1(k)+ Minute_1(k)/60+Second_1(k)/3600)/24,31)];
    end

    log1.time = date_1;
        
    meteo.time          = log1.time;
	meteo.temperature   = log1.meteo_temperature+273.15;
	meteo.rel_humidity  = log1.rel_humidity;
	meteo.pressure      = log1.pressure;
	meteo.wind_speed    = log1.wind_speed;
	meteo.wind_dirn     = log1.wind_dirn;
	meteo.visibility    = 0*ones(1,length(log1.meteo_temperature));
	meteo.global_rad    = 0*ones(1,length(log1.meteo_temperature));
	meteo.precipitation = log1.precipitation;
    	     
end