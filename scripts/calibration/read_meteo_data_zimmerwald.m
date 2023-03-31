function meteoData = read_meteo_data_zimmerwald(calibrationTool)
%==========================================================================
% NAME          | read_meteo_data_unibe.m
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014 // modified 01.2020
%               |
% ABSTRACT      | Function to read temperature and humidity values measured
%               | by the sensors on the top of ExWi building at Unibe. 
%               | If the local meteo station file does not exist,
%               | we try reading the anetz meteo data directly. 
%               |
% ARGUMENTS     | INPUTS: 1. calibrationTool:
%               |           - dateTime
%               |           - timeZone
%               |           - dateStr
%               |           - meteoFolder
%               |           - referenceTime
%               |           - meteoAnetzFolder
%               |           - anetzStnName
%               |
%               | OUTPUTS: 1. meteoData: structure containing the meteo
%               |            data read from exwi meteo station with the
%               |            folowing fields:
%               |           - dateNum
%               |           - dateTime
%               |           - air_temperature
%               |           - rel_humidity
%               |           - air_pressure
%               |           - precipitation
%               |           - tod
%               | 
% COMMENTS      | Time variations for this functions: TODO
%==========================================================================
%try
disp('Reading Zimmerwald Meteo')

dateStringMeteo=[calibrationTool.dateStr(3:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
meteoDataFile=[calibrationTool.meteoFolder 'L0' dateStringMeteo '.csv'];
%opts = detectImportOptions(meteoDataFile); % number of header lines which are to be ignored
%opts.DataLine = 3
%Tbl = dlmread(meteoDataFile, ' ' , 2, 0);
%disp(Tbl)


opts = detectImportOptions(meteoDataFile);%detectImportOptions(meteoDataFile); % number of header lines which are to be ignored
disp(opts)
opts.VariableNamesLine = 1 ; % row number which has variable names
%opts.ReadVariableNames = 'true'
%opts.VariableDescriptionsLine = 2
opts.DataLine = 3; % row number from which the actual data starts
%opts.VariableUnits = 2; % row number which has units specified
%disp('Reading Table')
opts.VariableNames(1) = {'time'};
opts.VariableNames(5) = {'air_TAmb_C'};
opts.VariableNames(11) = {'rel_humidity'};
opts.VariableNames(23) = {'pressure'};
opts.VariableNames(43) = {'wind_speed'};
opts.VariableNames(49) = {'wind_dirn'};
opts.VariableNames(69) = {'visibility'};
opts.VariableNames(79) = {'precipitation'};

%preview(meteoDataFile,opts);
T = readtable(meteoDataFile,opts);

%disp(meteoData(1).time)
%Transforming it into matlab structure

meteoData=table2struct(T);
disp('Reading Zimmerwald Meteo')
%disp(fieldnames(meteoData));

disp(meteoData(1).time)
disp('meteoData(5)')
disp(length(meteoData))

%meteoData(1).precipitation = meteoData(1).rain_accumulation;
for i = 1:length(meteoData)
    if meteoData(i).time.Year < 1900
        meteoData(i).time =  meteoData(i).time + years(2000);
    end
    %dateN = convertTo(meteoData(i).time, 'epochtime','Epoch','1970-01-01')/(24.*60.*60.);%dateN = datenum(meteoData(i).time);
    dateN = datenum(meteoData(i).time);
    meteoData(i).dateNum=dateN;
    meteoData(i).dateTime=meteoData(i).time;
    meteoData(i).temperature=meteoData(i).air_TAmb_C + 273.15;
    meteoData(i).air_temperature = meteoData(i).temperature;
    %disp(meteoData(i).air_temperature)
    
end
for i = 1:length(meteoData)
    meteoData(i).time = meteoData(i).dateNum;
    meteoData(i).air_pressure = meteoData(i).pressure;
end

end
