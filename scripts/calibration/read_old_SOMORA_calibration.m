function [data_FFT datalog_FFT] = read_old_SOMORA_calibration(mfile_name,rfile_name)

filename=mfile_name;

if exist(filename,'file')
      fid = fopen(filename,'r','ieee-le');
  if fid >= 3
    disp(['File ',filename,' loaded']);
    data_FFT = read_calibrated_FFT( fid );
    fclose(fid);
  elseif fid == -1
    disp(['Couldn''t open file ', filename]);
    data_FFT = [];
  end
else
    disp([ 'No SOMORA data for day ',num2str(date_to_invert)]);
    data_FFT = [];
end

%load calibrated data of log in datalog_FFT structure
% filename = [path.data,'r',num2str(date_to_invert),'.bin'];
filename=rfile_name;

if exist(filename,'file')
  fid = fopen(filename,'r','ieee-le');
  if fid >= 3
    disp(['File ',filename,' loaded']);
    %datalog = read_calibrated_FFT_log( fid );
    datalog_FFT = read_calibrated_FFT_log_v2( fid );
    fclose(fid);
  elseif fid == -1
    disp(['Couldn''t open file ', filename]);
    datalog_FFT = [];
  end


else
    disp([ 'No SOMORA data log for day ',num2str(date_to_invert)]);
    datalog_FFT = [];
end

function [record] = read_calibrated_FFT( fid )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading routine for SOMORA  FFT Version 1.0 Calibrated Data
% E.Maillard Barras (20100417)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        record = read_calibrated_FFT( fid )
%
% record: Matlab Structure Array
% fid:    Filehandle from fopen()
%
% Wichtig! Die Datei sollte unbedingt mit
% 
%        fid = fopen( 'name', 'r', 'ieee-le' )
%
% geöffnet werden, sonst wird das Intel-Binärformat der Datei nicht
% korrekt in das Binärformat des einlesenden Systems umgewandelt.
%

%fid = fopen(filename,'r','ieee-le');

for i=1:48
record(i).BT         = fread( fid, 16384, 'float32' ); % [K]
record(i).sigma_BT   = fread( fid, 16384, 'float32' ); % [K]
end

  
end


function [record] = read_calibrated_FFT_log_v2( fid )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading routine for SOMORA  FFT Version 1.0 Calibrated Data
% E.Maillard Barras (20100417)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        record = read_calibrated_FFT_log( fid )
%
% record: Matlab Structure Array
% fid:    Filehandle from fopen()
%
% Wichtig! Die Datei sollte unbedingt mit
% 
%        fid = fopen( 'name', 'r', 'ieee-le' )
%
% geöffnet werden, sonst wird das Intel-Binärformat der Datei nicht
% korrekt in das Binärformat des einlesenden Systems umgewandelt.
%

%fid = fopen(filename,'r','ieee-le');

for i=1:48 %not OK with 48 if t<48??
record(i).nbelement         = fread( fid, 1, 'int32' );
if isempty(record(i).nbelement)
    record(i).nbelement=0;
end
record(i).Tsys              = fread( fid, record(i).nbelement, 'float32' ); 
record(i).Twing             = fread( fid, record(i).nbelement, 'float32');
record(i).sigma_Tsys        = fread( fid, 1, 'float32' );
record(i).sigma_Twing       = fread( fid, 1, 'float32' );
record(i).Y                 = fread( fid, 1, 'int32' ); %9999
record(i).offset            = fread( fid, 1, 'int32'); % 0
record(i).angle_ant         = fread( fid, 1, 'float32');
record(i).sigma_angle_ant   = fread( fid, 1, 'float32');
record(i).date              = fread( fid, 1, 'int64'); 
record(i).t_mean            = fread( fid, 1, 'float32');
record(i).int_time          = fread( fid, 1, 'int32');
record(i).t_min             = fread( fid, 1, 'float32');
record(i).t_max             = fread( fid, 1, 'float32');
record(i).T_hot             = fread( fid, 1, 'float32');
record(i).T_room            = fread( fid, 1, 'float32');
record(i).T_out             = fread( fid, 1, 'float32');
record(i).T_win             = fread( fid, 1, 'float32');
record(i).T_res             = fread( fid, 1, 'bit32'); %9999
record(i).sigma_T_hot       = fread( fid, 1,'float32');
record(i).sigma_T_room      = fread( fid, 1,'float32');
record(i).sigma_T_out       = fread( fid, 1,'float32');
record(i).sigma_T_win       = fread( fid, 1,'float32');
record(i).sigma_T_res       = fread( fid, 1, 'bit32'); %9999
record(i).n_samples         = fread( fid, 1, 'int32' ); %length(data(t).t_avg)
record(i).cal_type          = fread( fid, 1, 'int32' ); %1=Planck

end

  
end
end
