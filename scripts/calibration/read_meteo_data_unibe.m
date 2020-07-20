function meteoData = read_meteo_data_unibe(calibrationTool)
%==========================================================================
% NAME          | get_meteo_data_unibe
% TYPE          | function
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014 // modified 01.2020
%               |
% ABSTRACT      | Function to read temperature and humidity values measured
%               | by the sensors on the top of ExWi building at Unibe
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
if calibrationTool.dateTime >= datetime(2017,08,10)
    % First reading the Meteo dataset for this day
    dateStringMeteo=[calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)];
    meteoDataFile=[calibrationTool.meteoFolder 'meteo_' dateStringMeteo '.csv'];
    
    % Transforming it into matlab structure
    T=readtable(meteoDataFile);
    meteoData=table2struct(T);
    
    meteoData(1).precipitation = meteoData(1).rain_accumulation;
    for i = 1:length(meteoData)
        dateN = datenum(meteoData(i).time,'yyyy-mm-ddTHH:MM:SS.FFFZ');
        meteoData(i).dateNum=dateN-datenum(1970,1,1);
        meteoData(i).dateTime=datetime(dateN,'ConvertFrom','datenum');
        %meteoData(i).dateTime=datetime(meteoData(i).time,'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSSSSSZ');
        meteoData(i).air_temperature=meteoData(i).air_temperature + calibrationTool.zeroDegInKelvin;
        meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
        % TODO Check units
        if i>1
            meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
        end
    end
elseif  (calibrationTool.dateTime > datetime(2017,01,01) && calibrationTool.dateTime < datetime(2017,08,10))
    disp('status of meteo data unkown between 01.01 and 09.08.2017') 
    meteoData = struct();
else
    dateStringMeteo=[calibrationTool.dateStr(3:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
    dateStringMeteoPrec=[calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
    baseName = [calibrationTool.meteoFolder calibrationTool.dateStr(1:4) '/' dateStringMeteo];
    meteoDataFileLog=[baseName '.log'];
    meteoDataFilePressure=[baseName 'pressure.log'];
    meteoDataFilePrecipitation=[calibrationTool.meteoFolder calibrationTool.dateStr(1:4) '/' dateStringMeteoPrec '_rainsensor.txt'];
    
    % Transforming it into matlab structure
    prec=readtable(meteoDataFilePrecipitation);
    precipitation=table2struct(prec);
    
    p=readtable(meteoDataFilePressure,'FileType','text');
    p.Properties.VariableNames = {'time' 'inst' 'mean' 'max','min'};
    pressure=table2struct(p);
%     %fmt = '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s  %*s  %*s  %*s  %*s  %*s  %*s  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f  %s  %s  %s  %*s %s  %s  %s  %s  %s  %s %*[^\n]';
%     fmt = '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f %*[^\n]';
%     %fmt = '%4d-%2d-%2d %2d:%2d %*s %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %f %f %f %f %s %f %f %f %f %f %f';
%     
%     %f  %f  %f  %f  %f  %f  %f  %f %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f';
% %     
% % %     % Transforming it into matlab structure
%      fileID = fopen(meteoDataFileLog,'r');
%      meteoFile = textscan(fileID,fmt,'Delimiter',' \b\t', 'EndOfLine','\n', 'TreatAsEmpty',{'/////' '//////'});
%      %'TreatAsEmpty',{'/////' '//////'}
% %     
%     for i = 1:144       
%         meteoData(i).dateTime=pressure(i).time;
%         meteoData(i).dateNum=meteoData(i).dateTime-datenum(1970,1,1);
%         %meteoData(i).dateTime=datetime(meteoData(i).time,'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSSSSSZ');
%         meteoData(i).air_temperature=meteoFile{1,9} + calibrationTool.zeroDegInKelvin;
%         meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
%         meteoData(i).precipitation = precipitation.Var2(i)
% %         % TODO Check units
% %         if i>1
% %             meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
% %         end
%     end
% % try
% 
%   if datenum(calibrationTool.Year, calibrationTool.Month, calibrationTool.Day) > 731616 & datenum(calibrationTool.Year, calibrationTool.Month, calibrationTool.Day) < 731673
%     % In dieser Periode hat der Windgeschwindigkeitssensor nicht funktioniert
% result = struct();
% [result.jahr ... % Datum und Zeit
%      result.monat ...
%      result.tag ...
%      result.stunde ...
%      result.minute ...
%      result.inst.temp ... % Temperatur
%      result.ave.temp ...
%      result.max.temp ...
%      result.min.temp ...
%      result.inst.relhum ... % Relative Feuchtigkeit
%      result.ave.relhum ...
%      result.max.relhum ...
%      result.min.relhum ...
%      result.inst.dewp ... % Taupunkt-Temperatur
%      result.ave.dewp ...
%      result.max.dewp ...
%      result.min.dewp ...
%      result.inst.airp ... % Luftdruck
%      result.ave.airp ...
%      result.max.airp ...
%      result.min.airp ...
%      result.inst.qnhp ... % QNH
%      result.ave.qnhp ...
%      result.max.qnhp ...
%      result.min.qnhp ...
%      result.inst.qffp ... % QFP
%      result.ave.qffp ...
%      result.max.qffp ...
%      result.min.qffp ...
%      result.inst.winddir2 ... % Windrichtung
%      result.ave.winddir2 ...
%      result.max.winddir2 ...
%      result.min.winddir2 ...
%      result.ave.winddir10 ... % Windrichtung
%      result.max.winddir10 ...
%      result.min.winddir10 ...
%      result.inst.sunrad ... % Sonnenstrahlung (Pyranometer)
%      result.ave.sunrad ...
%      result.max.sunrad ...
%      result.min.sunrad ...
%      result.inst.precip ... % Regen
%      result.sum.precip ... % Regen (summiert)
%      result.opvolt ... % Spannung Batterie
%      sta ... % Status FD12P
%      visibi1 ... % Optische Sichtweite
%      visibi2 ...
%      wecode1m ... % Wettercode instantan
%      wecode15m ... % Wettercode letzte 15min
%      wecode60m ... % Wettercode letzte 60min
%      rain1 ... % Regen (Sensor) instantan
%      rain2 ... % Regen (Sensor) summiert
%      snow ... % Schnee summiert
%     ]
%=textread(result.filename, '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s  %*s  %*s  %*s  %*s  %*s  %*s  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f  %s  %s  %s  %*s %s  %s  %s  %s  %s  %s',144);
%     result.inst.windsp2 = 999999*ones(144, 1); % Windgeschwindigkeit
%     result.ave.windsp2 = 999999*ones(144, 1);
%     result.max.windsp2 = 999999*ones(144, 1);
%     result.min.windsp2 = 999999*ones(144, 1);
%     result.ave.windsp10 = 999999*ones(144, 1); % Windgeschwindigkeit
%     result.max.windsp10 = 999999*ones(144, 1);
%     result.min.windsp10 = 999999*ones(144, 1);
%   else
%     [result.jahr ... % Datum und Zeit
%      result.monat ...
%      result.tag ...
%      result.stunde ...
%      result.minute ...
%      result.inst.temp ... % Temperatur
%      result.ave.temp ...
%      result.max.temp ...
%      result.min.temp ...
%      result.inst.relhum ... % Relative Feuchtigkeit
%      result.ave.relhum ...
%      result.max.relhum ...
%      result.min.relhum ...
%      result.inst.dewp ... % Taupunkt-Temperatur
%      result.ave.dewp ...
%      result.max.dewp ...
%      result.min.dewp ...
%      result.inst.airp ... % Luftdruck
%      result.ave.airp ...
%      result.max.airp ...
%      result.min.airp ...
%      result.inst.qnhp ... % QNH
%      result.ave.qnhp ...
%      result.max.qnhp ...
%      result.min.qnhp ...
%      result.inst.qffp ... % QFP
%      result.ave.qffp ...
%      result.max.qffp ...
%      result.min.qffp ...
%      result.inst.windsp2 ... % Windgeschwindigkeit
%      result.ave.windsp2 ...
%      result.max.windsp2 ...
%      result.min.windsp2 ...
%      result.ave.windsp10 ... % Windgeschwindigkeit
%      result.max.windsp10 ...
%      result.min.windsp10 ...
%      result.inst.winddir2 ... % Windrichtung
%      result.ave.winddir2 ...
%      result.max.winddir2 ...
%      result.min.winddir2 ...
%      result.ave.winddir10 ... % Windrichtung
%      result.max.winddir10 ...
%      result.min.winddir10 ...
%      result.inst.sunrad ... % Sonnenstrahlung (Pyranometer)
%      result.ave.sunrad ...
%      result.max.sunrad ...
%      result.min.sunrad ...
%      result.inst.precip ... % Regen
%      result.sum.precip ... % Regen (summiert)
%      result.opvolt ... % Spannung Batterie
%      sta ... % Status FD12P
%      visibi1 ... % Optische Sichtweite
%      visibi2 ...
%      wecode1m ... % Wettercode instantan
%      wecode15m ... % Wettercode letzte 15min
%      wecode60m ... % Wettercode letzte 60min
%      rain1 ... % Regen (Sensor) instantan
%      rain2 ... % Regen (Sensor) summiert
%      snow ... % Schnee summiert
%     ]=textread(result.filename, '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f  %s  %s  %s  %*s %s  %s  %s  %s  %s  %s',144);
%   end
%   
%   for i = 1:length(result.jahr)
%     if ~isempty(str2num(sta{i}))
%       result.sta(i,1) = str2num(sta{i});
%     else
%       result.sta(i,1) = 999;
%     end;
%     if ~isempty(str2num(visibi1{i}))
%       result.inst.visibi(i,1) = str2num(visibi1{i});
%     else
%       result.inst.visibi(i,1) = 99999;
%     end;
%     if ~isempty(str2num(visibi2{i}))
%       result.ave.visibi(i,1) = str2num(visibi2{i});
%     else
%       result.ave.visibi(i,1) = 99999;
%     end;
%     if ~isempty(str2num(wecode1m{i}))
%       result.wecode1m(i,1) = str2num(wecode1m{i});
%     else
%       result.wecode1m(i,1) = 999;
%     end;
%     if ~isempty(str2num(wecode15m{i}))
%       result.wecode15m(i,1) = str2num(wecode15m{i});
%     else
%       result.wecode15m(i,1) = 999;
%     end;
%     if ~isempty(str2num(wecode60m{i}))
%       result.wecode60m(i,1) = str2num(wecode60m{i});
%     else
%       result.wecode60m(i,1) = 999;
%     end;
%     if ~isempty(str2num(rain1{i}))
%       result.inst.rain(i,1) = str2num(rain1{i});
%     else
%       result.inst.rain(i,1) = 9999;
%     end;
%     if ~isempty(str2num(rain2{i}))
%       result.sum.rain(i,1) = str2num(rain2{i});
%     else
%       result.sum.rain(i,1) = 9999;
%     end;
%     if ~isempty(str2num(snow{i}))
%       result.sum.snow(i,1) = str2num(snow{i});
%     else
%       result.sum.snow(i,1) = 9999;
%     end;
%   end;
%   
%   result.tagdec = datenum(jahr, monat, tag, result.stunde, result.minute, 0); % Datum und Zeit im Matlab-Format
% 
%   result.leng = length(result.jahr); % L�nge der Vektoren
%   
% catch
% 
%   result.jahr = 999999;
%   result.monat = 999999;
%   result.tag = 999999;
%   result.stunde = 999999;
%   result.minute = 999999;
%   result.inst.temp = 999999;
%   result.ave.temp = 999999;
%   result.max.temp = 999999;
%   result.min.temp = 999999;
%   result.inst.relhum = 999999;
%   result.ave.relhum = 999999;
%   result.max.relhum = 999999;
%   result.min.relhum = 999999;
%   result.inst.dewp = 999999;
%   result.ave.dewp = 999999;
%   result.max.dewp = 999999;
%   result.min.dewp = 999999;
%   result.inst.airp = 999999;
%   result.ave.airp = 999999;
%   result.max.airp = 999999;
%   result.min.airp = 999999;
%   result.inst.qnhp = 999999;
%   result.ave.qnhp = 999999;
%   result.max.qnhp = 999999;
%   result.min.qnhp = 999999;
%   result.inst.qffp = 999999;
%   result.ave.qffp = 999999;
%   result.max.qffp = 999999;
%   result.min.qffp = 999999;
%   result.inst.windsp2 = 999999;
%   result.ave.windsp2 = 999999;
%   result.max.windsp2 = 999999;
%   result.min.windsp2 = 999999;
%   result.ave.windsp10 = 999999;
%   result.max.windsp10 = 999999;
%   result.min.windsp10 = 999999;
%   result.inst.winddir2 = 999999;
%   result.ave.winddir2 = 999999;
%   result.max.winddir2 = 999999;
%   result.min.winddir2 = 999999;
%   result.ave.winddir10 = 999999;
%   result.max.winddir10 = 999999;
%   result.min.winddir10 = 999999;
%   result.inst.sunrad = 999999;
%   result.ave.sunrad = 999999;
%   result.max.sunrad = 999999;
%   result.min.sunrad = 999999;
%   result.inst.precip = 999999;
%   result.sum.precip = 999999;
%   result.opvolt = 999999;
%   result.sta = 999999;
%   result.inst.visibi = 999999;
%   result.ave.visibi = 999999;
%   result.wecode1m = 999999;
%   result.wecode15m = 999999;
%   result.wecode60m = 999999;
%   result.inst.rain = 999999;
%   result.sum.rain = 999999;
%   result.sum.snow = 999999;
% 
%   result.tagdec = 9999999;
% 
%   result.leng = 0;
%     
% end;


%     
    disp('format of meteo data changed before the 10th of August 2017') 
    meteoData = struct();
end
