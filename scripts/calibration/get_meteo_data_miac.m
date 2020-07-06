function meteo = get_meteo_data_miac(logFile)
%==========================================================================
% NAME          | get_meteo_data_miac
% TYPE          | function
% AUTHOR(S)     | Adapted from 'insert_meteo_data_into_mysql_v3', by F. Schranz
% CREATION      | 2020-06-12
%               |
% ABSTRACT      | Function to read meteo data from the MIAWARA-C weather
%               | station when it is saved in the MIAWARA-C log file.
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

if ~isempty(logFile) 
	meteo.time          = logFile.time;
	meteo.temperature   = logFile.temperature+273.15;
	meteo.rel_humidity  = logFile.rel_humidity;
	meteo.pressure      = logFile.pressure;
	meteo.wind_speed    = logFile.wind_speed;
	meteo.wind_dirn     = logFile.wind_dirn;
	meteo.visibility    = 0*ones(1,length(logFile.temperature));
	meteo.global_rad    = 0*ones(1,length(logFile.temperature));
	meteo.precipitation = logFile.precipitation;
	
else 
	meteo.time          = [];
	meteo.temperature   = [];
	meteo.rel_humidity  = [];
	meteo.pressure      = [];
	meteo.wind_speed    = [];
	meteo.wind_dirn     = [];
	meteo.visibility    = [];
	meteo.global_rad    = [];
	meteo.precipitation = [];
	
end