function read_gromos_file( year, month );

% checking if the year and month are valid for GROMOS series

if (year < 1994) | (year > 2012) | (month < 1) | (month > 12)
  error('Invalid year or month!');
end;

% go to folder where data-files are

folder = ['/storage/atmosphere/instruments/gromos/level1/spectra/o3', num2str( year )];
%eval( sprintf('cd %s', folder) );
data = struct();
meteoData= struct();
% getting month and date as stored in the names of GROMOS files

switch month
  case 1
    mon = 'jan';
    length = 31;
  case 2
    mon = 'feb';
    length = 28;
  case 3
    mon = 'mar';
    length = 31;
  case 4
    mon = 'apr';
    length = 30;
  case 5
    mon = 'may';
    length = 31;
  case 6
    mon = 'jun';
    length = 30;
  case 7
    mon = 'jul';
    length = 31;
  case 8
    mon = 'aug';
    length = 31;
  case 9
    mon = 'sep';
    length = 30;
  case 10
    mon = 'oct';
    length = 31;
  case 11
    mon = 'nov';
    length = 30;
  case 12
    mon = 'dec';
    length = 31;
end;

if ((floor(year/4) - year/4) == 0 ) && (month == 2)
  length = 29;
end;

y = mod( year, 100 );
if y >= 10
  ye = num2str( y );
else
  ye = strcat( '0', num2str( y ) );
end;


% if month < 10
%   file_to_write = strcat('/scratch/GROSOM/Level1/GROMOS/FB/gromos_', num2str(year), '_0', num2str(month), '.txt');
% else
%   file_to_write = strcat('/scratch/GROSOM/Level1/GROMOS/FB/gromos_', num2str(year), '_', num2str(month), '.txt');
% end;
% fidw = fopen( file_to_write, 'w' );

%clear zeile;

% reading of all monthly files and storage into the variable 'zeile'

for i = 1: length

% obtaining the name of the file to be read, and opening of that file

  if i >= 10
    day = num2str( i );
  else
    day = strcat( '0', num2str( i ) );
  end;
  file = [folder '/' strcat('o3', mon , ye, '.d', day)];
  if exist( file )
    fid = fopen( file );
  else
    fid = fopen( upper(file) );
  end;

% if the file is not found, it should go to the next one

  if fid < 0
   disp ( ['file ' file ' not found']);
   continue;
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
      %process = 'yes';
      head = sscanf(row, '%s');
      x = sscanf(row(44:76), '%f\t%f\t%f\t%f\t%f\t%f\t%f\t');
      data(entries).date             = head(9:17);
      data(entries).t_first          = head(18:25);
      data(entries).t_last           = head(1:8);
      %time_first       = convert_to_mysql_date( date, t_first );
      %time_last        = convert_to_mysql_date( date, t_last );
      data(entries).nt               = x(2);
      data(entries).nt_var           = x(3);
      data(entries).t_hot            = x(4);
      data(entries).t_cold           = x(5);
      if channels_ready
        if (mean(channel) < 5) | (data(entries).nt > 5000) | (data(entries).nt < 1000)
	        process = 'no';
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
  disp( sprintf('day %d of month %d in %d is O.K.\n', i, month, year) );
end

fclose( 'all' );
