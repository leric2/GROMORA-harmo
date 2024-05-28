function [data, meteoData, calibrationTool] = read_level1_FB_GROMORA(calibrationTool)
%==========================================================================
% NAME          | read_level1_FB_GROMORA.m
% TYPE          | function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2023
%               |
% ABSTRACT      | Function to read a level1 calibrated file from the old FB
%               | data. 
%               |
% ARGUMENTS     | INPUTS:   
%               |  
%               | OUTPUTS:
%               |
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the different dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = struct();
meteoData= struct();

% checking if the Year and month are valid for GROMOS series
% 
if (calibrationTool.Year < 1989) | (calibrationTool.Year > 2012) | (calibrationTool.Month < 1) | (calibrationTool.Month > 12)
  error('Invalid year or month!');
end

% go to folder where data-files are
if (calibrationTool.Year < 1994)
    folder = '/storage/atmosphere/instruments/gromos/level1/spectra/old_spectra';
else
    folder = ['/storage/atmosphere/instruments/gromos/level1/spectra/o3', num2str( calibrationTool.Year )];
end
%eval( sprintf('cd %s', folder) );
data = struct();

% getting month and date as stored in the names of GROMOS files

switch calibrationTool.Month
  case 1
    mon = 'jan';
    month_length = 31;
  case 2
    mon = 'feb';
    month_length = 28;
  case 3
    mon = 'mar';
    month_length = 31;
  case 4
    mon = 'apr';
    month_length = 30;
  case 5
    mon = 'may';
    month_length = 31;
  case 6
    mon = 'jun';
    month_length = 30;
  case 7
    mon = 'jul';
    month_length = 31;
  case 8
    mon = 'aug';
    month_length = 31;
  case 9
    mon = 'sep';
    month_length = 30;
  case 10
    mon = 'oct';
    month_length = 31;
  case 11
    mon = 'nov';
    month_length = 30;
  case 12
    mon = 'dec';
    month_length = 31;
end

if ((floor(calibrationTool.Year/4) - calibrationTool.Year/4) == 0 ) && (calibrationTool.Month == 2)
  month_length = 29;
end

y = mod( calibrationTool.Year, 100 );
if y >= 10
  ye = num2str( y );
else
  ye = strcat( '0', num2str( y ) );
end


% if calibrationTool.Month < 10
%   file_to_write = strcat('/scratch/GROSOM/Level1/GROMOS/FB/gromos_', num2str(calibrationTool.Year), '_0', num2str(calibrationTool.Month), '.txt');
% else
%   file_to_write = strcat('/scratch/GROSOM/Level1/GROMOS/FB/gromos_', num2str(calibrationTool.Year), '_', num2str(calibrationTool.Month), '.txt');
% end;
% fidw = fopen( file_to_write, 'w' );

%clear zeile;

% reading of all monthly files and storage into the variable 'zeile'

%for i = 1: month_length
i = calibrationTool.Day;
% obtaining the name of the file to be read, and opening of that file

  if i >= 10
    day = num2str( i );
  else
    day = strcat( '0', num2str( i ) );
  end
  file = [folder '/' strcat('o3', mon , ye, '.d', day)];
  if exist( file )
    fid = fopen( file );
  else
    file = [folder '/' upper(strcat('o3', mon , ye, '.d', day))];
    fid = fopen( file );
  end

% filename:
calibrationTool.filenameLevel1a = file;


% if the file is not found, it should go to the next one
  if fid < 0
   error( ['file ' file ' not found']);
  end

% reading of every line from the file; if the line contains channels, the values
% are stored in the variable 'channel'. When it comes to the header, the haeader
% data are stored at the begining of the variable 'zeile'
  entries = 1;
  j = 1;
  channels_ready = 0;

  row = fgetl( fid );
  while not( isempty(row) ) && not( size(row, 2) < 64 ) && not( feof( fid) )
    if (row(1) == ' ') && (size(row, 2) == 64) % when the line begins with ' ' and contains numbers, these are channels
        channel( (j-1)*8+1:(j-1)*8+8 ) = sscanf( row, '%f' )';
        data(entries).Tb = channel;
        j = j + 1;
	if size(channel, 2) == 48
	  channels_ready = 1;
    end
    elseif (row(3) == ':') && (row(6) == ':') % when the line begins with time in format 'xx:yy:zz', this is header
      data(entries).process = 'yes';
      head = sscanf(row, '%s');
      x = sscanf(row(44:76), '%f\t%f\t%f\t%f\t%f\t%f\t%f\t');
      data(entries).day             = str2num(head(9:10));
      data(entries).datestr       = [head(9:12) lower(head(13:14)) head(15:17)];
      data(entries).month =  calibrationTool.Month; 
      data(entries).year             = calibrationTool.Year;
      %data(entries).datenum = datenum(data(entries).day , ,);
      %data(entries).dateTime = datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day, 'TimeZone',calibrationTool.timeZone);
      data(entries).beginningDay = datenum(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day) - calibrationTool.referenceTime;
      
      if entries == 1 && str2num(head(18:19)) > 0
            data(entries).first_sky_time = datenum(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day-1,str2num(head(18:19)),str2num(head(21:22)),str2num(head(24:25))) - calibrationTool.referenceTime;
      else
            data(entries).first_sky_time = datenum(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day,str2num(head(18:19)),str2num(head(21:22)),str2num(head(24:25))) - calibrationTool.referenceTime;
      end
      data(entries).time_min = datenum(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day,str2num(head(18:19)),str2num(head(21:22)),str2num(head(24:25))) - calibrationTool.referenceTime;
      data(entries).last_sky_time = datenum(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day,str2num(head(1:2)),str2num(head(4:5)),str2num(head(7:8))) - calibrationTool.referenceTime;
      data(entries).mean_sky_time =(data(entries).first_sky_time+ data(entries).last_sky_time)/2;
      data(entries).dateTime =datetime(data(entries).mean_sky_time+calibrationTool.referenceTime, 'ConvertFrom','datenum', 'TimeZone',calibrationTool.timeZone);
      data(entries).t_first          = head(18:25);
      data(entries).t_last           = head(1:8);
      data(entries).calibration_time_duration = datetime(data(entries).last_sky_time+calibrationTool.referenceTime, 'ConvertFrom','datenum', 'TimeZone',calibrationTool.timeZone)-datetime(data(entries).first_sky_time+calibrationTool.referenceTime, 'ConvertFrom','datenum', 'TimeZone',calibrationTool.timeZone);
      data(entries).calibration_time = seconds(data(entries).calibration_time_duration);
      data(entries).some_time_calibration = head(26:33);
      data(entries).a_flag           = str2num(head(34)) ;
      data(entries).another_flag           = str2num(head(35));
      %time_first       = convert_to_mysql_date( date, t_first );
      %time_last        = convert_to_mysql_date( date, t_last );
      data(entries).TNoise               = x(2);
      data(entries).stdTNoise           = x(3);
      data(entries).THot            = x(4);
      data(entries).Tcold           = x(5);
      data(entries).mean_sky_elevation_angle  = x(6);
      data(entries).opacity           = x(7);
      data(entries).flags           = 1*[str2num(head(34))==0 str2num(head(35))==0];
      if channels_ready
        if (mean(data(entries).Tb) < 5) | (data(entries).TNoise > 5000) | (data(entries).TNoise < 1000)
	         data(entries).process = 'no';
    	end
        %channels_string  = sprintf('%9.2f\t', channel(1:48));
        %zeile = sprintf( '%s\t%d\t%s\t%s\t%9.2f\t%8.2f\t%6.1f\t%6.1f\t%s\t%s\n', '\N', 1, t_first, t_last, nt, nt_var, t_hot, t_cold, process, channels_string);
        %fwrite( fidw, zeile ); % writting the line in the file
        j = 1;
	    channels_ready = 0;
        entries=entries+1;
      end
    else
      disp(['Unknown format of line: ', row]);
      row = fgetl( fid );
      continue;
    end
    row = fgetl( fid );
  end
  fclose(fid);
  fprintf('day %d of month %d in %d is O.K.\n', i, calibrationTool.Month, calibrationTool.Year);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding additional data for successful integration within GROMORA standard
% routine
%data.dateTime = data.t_first
%data.meanTime = 
for i =1:length(data)
    data(i).stdTb = -9999*ones(length(data(i).Tb));
    data(i).noiseLevel = -9999;
    data(i).frequencies = calibrationTool.LOFreqTot + calibrationTool.intermediate_freq;
    data(i).intermediate_freq = calibrationTool.intermediate_freq;
    data(i).TWindow = 293;
    data(i).TRoom =-9999;
    data(i).VGunn =-9999;
    data(i).time_of_day = (data(i).mean_sky_time - data(i).beginningDay);
    data(i).number_of_hot_spectra = data(i).calibration_time/(3*calibrationTool.cycleDurationHot);  %[];
    data(i).number_of_cold_spectra = data(i).calibration_time/(3*calibrationTool.cycleDurationCold);  %[];
    data(i).number_of_sky_spectra = data(i).calibration_time/(3*calibrationTool.cycleDurationSky);  %[];
    if isempty(data(i).dateTime)
        data(i) = [];
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding meteo data
meteo = read_meteo_data_unibe(calibrationTool);
if isempty(fieldnames(meteo))
    meteoData = struct();
    meteoData.dateNum = nan;
else
meteoData = struct();
meteoData.dateNum = [meteo.dateNum];
meteoData.dateTime = [meteo.dateTime];
meteoData.air_pressure = [meteo.air_pressure];
meteoData.air_temperature = [meteo.air_temperature];
meteoData.relative_humidity = [meteo.rel_humidity];
meteoData.precipitation = [meteo.precipitation];
end
%     if strcmp(gName,'meteo')
%         vNames = {ncinfo(filename,gName).Variables.Name};
%         %dimNames = {ncinfo(filename,gName).Dimensions.Name};
%         %         size=zeros(1,length(dimNames));
%         %         for i=1:length(dimNames)
%         %             size(i) = ncinfo(filename,fullfile(gName,dimNames{i})).Size;
%         %         end
%         
%         for v=1:length(vNames)
%             varName = vNames{v};
%             
%             if length(ncinfo(filename,fullfile(gName,varName)).Dimensions)==1
%                 if strcmp(ncinfo(filename,fullfile(gName,varName)).Dimensions.Name,'time')
%                     meteoData.(varName) = ncread(filename,fullfile(gName,varName))';
%                 else
%                     meteoData(t).(varName) = ncread(filename,fullfile(gName,varName),1,Inf);
%                 end
%             else
%                 meteoData(t).(varName)= ncread(filename,fullfile(gName,varName),[1,t],[Inf,1])';
%             end
%         end
%     else
%         vNames = {ncinfo(filename,gName).Variables.Name};
%         dimNames = {ncinfo(filename,gName).Dimensions.Name};
%         size=zeros(1,length(dimNames));
%         for i=1:length(dimNames)
%             size(i) = ncinfo(filename,fullfile(gName,dimNames{i})).Size;
%         end
%         
%         % Intialize our structure with time
%         for t = 1:size(1)
%             varName = vNames{1};
%             data(t).(varName) = ncread(filename,fullfile(gName,varName),t,1);
%         end
%         dateT = num2cell(datetime([data.time] + calibrationTool.referenceTime,'ConvertFrom','datenum','TimeZone',calibrationTool.timeZone));
%         [data.dateTime] = dateT{:};
%         %if strcmp('time',dimName)
%         for v=2:length(vNames)
%             varName = vNames{v};
%             
%             if length(ncinfo(filename,fullfile(gName,varName)).Dimensions)==1
%                 if strcmp(ncinfo(filename,fullfile(gName,varName)).Dimensions.Name,'time')
%                     dat=num2cell(ncread(filename,fullfile(gName,varName)));
%                     [data.(varName)] = dat{:};
%                 else
%                     for t = 1:size(1)
%                         data(t).(varName) = ncread(filename,fullfile(gName,varName),1,Inf)';
%                     end
%                 end
%             else
%                 %dat=num2cell(ncread(filename,fullfile(gName,varName)))';
%                 for t = 1:size(1)
%                     data(t).(varName)= ncread(filename,fullfile(gName,varName),[1,t],[Inf,1])';
%                 end
%             end
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding some attributes
calibrationTool.logFile.raw_file_warning='';
calibrationTool.logFile.comment='';
calibrationTool.logFile.raw_file_comment='';
calibrationTool.logFile.raw_data_software_version='';
calibrationTool.logFile.calibration_version='';
calibrationTool.logFile.creation_date_level1a='';
calibrationTool.logFile.raw_data_software_version='';
calibrationTool.logFile.filenameLevel1a = file;
calibrationTool.flagVectorLength = 2;
calibrationTool.logFile.rawFilename = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding TC data
if calibrationTool.doTippingCurve
    file_TC = [folder '/' strcat('tip', mon , ye, '.d', day)];
  if exist( file_TC )
    fid = fopen(file_TC) ;
  else
    fid = fopen(upper(file_TC));
  end

% if the file is not found, it should go to the next one
  if fid < 0
   error( ['file ' file_TC ' not found']);
  end

  entriesTC = 1;
  row = fgetl( fid );
  while not( isempty(row) ) && not( size(row, 2) < 64 ) && not( feof( fid) )
      head = sscanf(row, '%s');
      calibrationTool.logFile.TC(entriesTC).Tb   = sscanf(row(22:end), '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t');
      calibrationTool.logFile.TC(entriesTC).day             = str2num(head(9:10));
      calibrationTool.logFile.TC(entriesTC).datestr       = [head(9:12) lower(head(13:14)) head(15:17)];
      calibrationTool.logFile.TC(entriesTC).month =  calibrationTool.Month; 
      calibrationTool.logFile.TC(entriesTC).year             = calibrationTool.Year;
      %TC(entriesTC).datenum = datenum(TC(entriesTC).day , ,);
      %TC(entriesTC).dateTime = datetime(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day, 'TimeZone',calibrationTool.timeZone);
     
      calibrationTool.logFile.TC(entriesTC).time = datenum(calibrationTool.Year,calibrationTool.Month,calibrationTool.Day,str2num(head(1:2)),str2num(head(4:5)),str2num(head(7:8))) - calibrationTool.referenceTime;
      calibrationTool.logFile.TC(entriesTC).dateTime =datetime(calibrationTool.logFile.TC(entriesTC).time+calibrationTool.referenceTime, 'ConvertFrom','datenum', 'TimeZone',calibrationTool.timeZone);
      calibrationTool.logFile.TC(entriesTC).frequency = reshape(-9999*ones(length(calibrationTool.logFile.TC(entriesTC).Tb)),1,[]);
      calibrationTool.logFile.TC(entriesTC).THot_calib = mean([data.THot]);
      %calibrationTool.logFile.TC(j).tau_tc = TC.tauCalib(j);
      calibrationTool.logFile.TC(entriesTC).tipping_angle = [45 41.25 37.71 34.33 31.09 27.95 24.89 21.92 19];%reshape(-9999*ones(length(calibrationTool.logFile.TC(entriesTC).Tb)),1,[]);
        if isfield(meteoData,'air_temperature')
            Teff = nanmean([meteoData.air_temperature])-calibrationTool.TC.deltaT;
        else
            %disp('we said, no meteo data found so lets make a guess for Tair (10 degC)');
            Teff = 283 - calibrationTool.TC.deltaT;
        end
            calibrationTool.logFile.TC(entriesTC).am = 1./sind(calibrationTool.logFile.TC(entriesTC).tipping_angle);
            tau = calibrationTool.TC.tauInitTC;
            it_max = calibrationTool.TC.maxIterTC;
            N_it = 0;
            off = 10;
            calibrationTool.logFile.TC(entriesTC).converged = 1;
            while abs(off) > calibrationTool.TC.offsetTC
                % Computation of T_cold with this tau
                TCold=calibrationTool.backgroundMWTb.*exp(-tau)+Teff.*(1-exp(-tau));
                %disp(T_cold)
                % Computation of Tb_theta for each angle with this T_cold
                Tb_theta=calibrationTool.logFile.TC(entriesTC).Tb;
                
                % For every angle, compute the corresponding tau_theta
                tau_theta=log((Teff-calibrationTool.backgroundMWTb)./(Teff-Tb_theta));
                %disp(tau_theta)
                
                % Linear fit for tau_theta vs airmass factor
                fit=polyfit( calibrationTool.logFile.TC(entriesTC).am,tau_theta,1);
                
                % The new tau is the slope of the linear regression
                off = fit(2);
                tau = fit(1);
                N_it=N_it+1;
                
                if N_it>it_max
                    %tau=NaN;
                    %disp('Failed to converge')
                    calibrationTool.logFile.TC(entriesTC).converged = 0;
                    break
                end
            end
            
            %calibrationTool.logFile.TC(j).Tb_fromTCLoads = calibrationTool.TCold + (calibrationTool.logFile.TC(j).THot - calibrationTool.TCold) .* (calibrationTool.logFile.TC(j).sky - calibrationTool.logFile.TC(j).cold)./(calibrationTool.logFile.TC(j).hot - calibrationTool.logFile.TC(j).cold);
            %tau_slant = log((Teff-calibrationTool.backgroundMWTb)./(Teff-calibrationTool.logFile.TC(j).Tb_fromTCLoads));
            
            calibrationTool.logFile.TC(entriesTC).tau_tc = real(tau);
       
      row = fgetl( fid );
      entriesTC = entriesTC+1;
  end
    
end

% if sublevel == 1
% calibrationTool.logFile.rawFilename=ncreadatt(filename,'/','raw_filename');
% calibrationTool.logFile.rawData=ncreadatt(filename,'/','raw_data');
% end
% 
% % Coordinate variables, directly adding the attributes
% calibrationTool.timeUnit = ncreadatt(filename,'/spectrometer1/time','units');
% calibrationTool.timeCalendar= ncreadatt(filename,'/spectrometer1/time','calendar');
% %calibrationTool.timeZone = ncreadatt(filename,'/spectrometer1/time','time_zone');
% 
% meteoData.dateTime = datetime(meteoData.time + calibrationTool.referenceTime,'ConvertFrom','datenum');
% meteoData.dateTime.TimeZone = calibrationTool.timeZone;


disp(['File read : ' file])
end


