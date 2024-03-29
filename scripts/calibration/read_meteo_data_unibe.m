function meteoData = read_meteo_data_unibe(calibrationTool)
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
try
    if calibrationTool.dateTime >= datetime(2017,08,10, 'TimeZone', calibrationTool.timeZone)
        % First reading the Meteo dataset for this day
        dateStringMeteo=[calibrationTool.dateStr(1:4) '-' calibrationTool.dateStr(6:7) '-' calibrationTool.dateStr(9:10)];
        meteoDataFile=[calibrationTool.meteoFolder 'meteo_' dateStringMeteo '.csv'];
        
        % Transforming it into matlab structure
        T=readtable(meteoDataFile);
        meteoData=table2struct(T);
        
        meteoData(1).precipitation = meteoData(1).rain_accumulation;
        for i = 1:length(meteoData)
            dateN = datenum(meteoData(i).time,'yyyy-mm-ddTHH:MM:SS.FFFZ');
            meteoData(i).dateNum=dateN-calibrationTool.referenceTime;
            meteoData(i).dateTime=datetime(dateN,'ConvertFrom','datenum');
            meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
            %meteoData(i).dateTime=datetime(meteoData(i).time,'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSSSSSZ');
            meteoData(i).air_temperature=meteoData(i).air_temperature + calibrationTool.zeroDegInKelvin;
            meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
            % TODO Check units
            if i>1
                meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
            end
        end
%     elseif  (calibrationTool.dateTime > datetime(2017,01,01, 'TimeZone', calibrationTool.timeZone) && calibrationTool.dateTime < datetime(2017,08,10, 'TimeZone', calibrationTool.timeZone))
%         disp('status of meteo data unkown between 01.01 and 09.08.2017')
%         meteoData = struct();
    elseif calibrationTool.dateTime < datetime(2017,08,10, 'TimeZone', calibrationTool.timeZone) && calibrationTool.dateTime >= datetime(2010,05,12, 'TimeZone', calibrationTool.timeZone)
        dateStringMeteo=[calibrationTool.dateStr(3:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        dateStringMeteoPrec=[calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        baseName = [calibrationTool.meteoFolder dateStringMeteo];
        meteoDataFileLog=[baseName '.log'];
        meteoDataFilePressure=[baseName 'pressure.log'];
        meteoDataFilePrecipitation=[calibrationTool.meteoFolder dateStringMeteoPrec '_rainsensor.txt'];
        
        % Transforming it into matlab structure
        precipitation=readtable(meteoDataFilePrecipitation,'FileType','text','TreatAsEmpty',{'//////'});
        precipitation.Properties.VariableNames = {'dateTime','p1','p2'};
        precipitation.dateTime.TimeZone = calibrationTool.timeZone;
        %precipitation=table2struct(prec);
        
        pressure=readtable(meteoDataFilePressure,'FileType','text','TreatAsEmpty',{'//','///','////','/////','//////','//////'});
        pressure.Properties.VariableNames = {'dateTime' 'inst' 'mean' 'max','min'};
        pressure.dateTime.TimeZone = calibrationTool.timeZone;
        %pressure=table2struct(p);
        %     %fmt = '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s  %*s  %*s  %*s  %*s  %*s  %*s  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f  %s  %s  %s  %*s %s  %s  %s  %s  %s  %s %*[^\n]';
        %fmt = '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f %*[^\n]';
        %     %fmt = '%4d-%2d-%2d %2d:%2d %*s %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %f %f %f %f %s %f %f %f %f %f %f';
        %
        %     %f  %f  %f  %f  %f  %f  %f  %f %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f';
        % %
        % % %     % Transforming it into matlab structure
        
        meteoFile=readtable(meteoDataFileLog,'FileType','text','TreatAsEmpty',{'//','///','////','/////','//////'});
        %fileID = fopen(meteoDataFileLog,'r');
        %meteoFile.year = textscan(fileID,fmt,'Delimiter',' \b\t', 'EndOfLine','\n', 'TreatAsEmpty',{'/////' '//////'});
        %      %'TreatAsEmpty',{'/////' '//////'}
        % %
        
        dt = minutes(10);
        meteoData=struct();
        for i = 1:height(meteoFile)
            meteoRow = meteoFile(i,:);
            meteoData(i).dateTime=meteoRow.Var1;
            meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
            meteoData(i).dateNum=datenum(meteoData(i).dateTime)-calibrationTool.referenceTime;
            meteoData(i).air_temperature=meteoRow.Var4 + calibrationTool.zeroDegInKelvin;
            meteoData(i).tod = 24*(meteoData(i).dateTime-meteoData(1).dateTime);
            
            %rowPrec = find(isbetween(precipitation.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt));
            %meteoData(i).precipitation = sum(precipitation(rowPrec,:).p1,'omitnan');
            meteoData(i).precipitation = nanmean(precipitation.p1(isbetween(precipitation.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt),:));;
            meteoData(i).rel_humidity=meteoRow.Var8;
            
            meteoData(i).air_pressure=nanmean(pressure(isbetween(pressure.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt),:).mean);
            %         % TODO Check units
            %         if i>1
            %             meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
            %         end
        end
       elseif calibrationTool.dateTime < datetime(2010,05,12, 'TimeZone', calibrationTool.timeZone) && calibrationTool.dateTime >= datetime(2008,05,22, 'TimeZone', calibrationTool.timeZone)
        dateStringMeteo=[calibrationTool.dateStr(3:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        dateStringMeteoPrec=[calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        baseName = [calibrationTool.meteoFolder dateStringMeteo];
        meteoDataFileLog=[baseName '.log'];
        meteoDataFilePrecipitation=[calibrationTool.meteoFolder dateStringMeteoPrec '_rainsensor.txt'];
        
        % Transforming it into matlab structure
        precipitation=readtable(meteoDataFilePrecipitation,'FileType','text','TreatAsEmpty',{'//////'});
        precipitation.Properties.VariableNames = {'dateTime','p1','p2'};
        precipitation.dateTime.TimeZone = calibrationTool.timeZone;
        %precipitation=table2struct(prec);
        
        %pressure=readtable(meteoDataFilePressure,'FileType','text','TreatAsEmpty',{'//','///','////','/////','//////','//////'});
        %pressure.Properties.VariableNames = {'dateTime' 'inst' 'mean' 'max','min'};
        %pressure.dateTime.TimeZone = calibrationTool.timeZone;
        %pressure=table2struct(p);
        %     %fmt = '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s  %*s  %*s  %*s  %*s  %*s  %*s  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f  %s  %s  %s  %*s %s  %s  %s  %s  %s  %s %*[^\n]';
        %fmt = '%4d-%2d-%2d %2d:%2d %*s %f  %f  %f  %f  %f  %f  %f  %f  %f %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f %*[^\n]';
        %     %fmt = '%4d-%2d-%2d %2d:%2d %*s %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %f %f %f %f %s %f %f %f %f %f %f';
        %
        %     %f  %f  %f  %f  %f  %f  %f  %f %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f';
        % %
        % % %     % Transforming it into matlab structure
        
        meteoFile=readtable(meteoDataFileLog,'FileType','text','TreatAsEmpty',{'//','///','////','/////','//////'});
        %fileID = fopen(meteoDataFileLog,'r');
        %meteoFile.year = textscan(fileID,fmt,'Delimiter',' \b\t', 'EndOfLine','\n', 'TreatAsEmpty',{'/////' '//////'});
        %      %'TreatAsEmpty',{'/////' '//////'}
        % %
        
        dt = minutes(10);
        meteoData=struct();
        for i = 1:height(meteoFile)
            meteoRow = meteoFile(i,:);
            meteoData(i).dateTime=meteoRow.Var1;
            meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
            meteoData(i).dateNum=datenum(meteoData(i).dateTime)-calibrationTool.referenceTime;
            meteoData(i).air_temperature=meteoRow.Var4 + calibrationTool.zeroDegInKelvin;
            meteoData(i).tod = hours(meteoData(i).dateTime-meteoData(1).dateTime);
            
            %rowPrec = find(isbetween(precipitation.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt));
            %meteoData(i).precipitation = sum(precipitation(rowPrec,:).p1,'omitnan');
            meteoData(i).precipitation = nanmean(precipitation.p1(isbetween(precipitation.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt),:));
            meteoData(i).rel_humidity=meteoRow.Var8;
            
            meteoData(i).air_pressure= nan;
            %         % TODO Check units
            %         if i>1
            %             meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
            %         end
        end
    elseif calibrationTool.dateTime < datetime(2008,05,22, 'TimeZone', calibrationTool.timeZone) && calibrationTool.dateTime >= datetime(1997,06,20, 'TimeZone', calibrationTool.timeZone)
        dateStringMeteo=[calibrationTool.dateStr(3:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        baseName = [calibrationTool.meteoFolder dateStringMeteo];
        meteoDataFileLog=[baseName '.log'];
        
        meteoFile=readtable(meteoDataFileLog,'FileType','text','TreatAsEmpty',{'//','///','////','/////','//////'});

        
        dt = minutes(10);
        meteoData=struct();
        for i = 1:height(meteoFile)
            meteoRow = meteoFile(i,:);
            meteoData(i).dateTime=meteoRow.Var1;
            meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
            meteoData(i).dateNum=datenum(meteoData(i).dateTime)-calibrationTool.referenceTime;
            meteoData(i).air_temperature=meteoRow.Var4 + calibrationTool.zeroDegInKelvin;
            meteoData(i).tod = hours(meteoData(i).dateTime-meteoData(1).dateTime);
            
            %rowPrec = find(isbetween(precipitation.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt));
            %meteoData(i).precipitation = sum(precipitation(rowPrec,:).p1,'omitnan');
            meteoData(i).precipitation = nan;
            meteoData(i).rel_humidity=meteoRow.Var8;
            
            meteoData(i).air_pressure= nan;
            %         % TODO Check units
            %         if i>1
            %             meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
            %         end
        end
    elseif calibrationTool.dateTime < datetime(1997,06,20, 'TimeZone', calibrationTool.timeZone)
        disp('No Exwi Meteo data, using Bollwerk station.')
        dateStringMeteo=['bollwerk_' calibrationTool.dateStr(1:4) calibrationTool.dateStr(6:7) calibrationTool.dateStr(9:10)];
        baseName = [calibrationTool.meteoFolderBollwerk dateStringMeteo];
        meteoDataFileLog=[baseName '.txt'];
        
        meteoFile=readtable(meteoDataFileLog,'FileType','text','TreatAsEmpty',{'-9999','//','///','////','/////','//////'}, 'HeaderLines', 22);

        dt = minutes(10);
        meteoData=struct();
        for i = 1:height(meteoFile)
            meteoRow = meteoFile(i,:);
            meteoData(i).dateTime=datetime(meteoRow.Var1,meteoRow.Var2, meteoRow.Var3, meteoRow.Var4, meteoRow.Var5, 0);
            meteoData(i).dateTime.TimeZone = calibrationTool.timeZone;
            meteoData(i).dateNum=datenum(meteoData(i).dateTime)-calibrationTool.referenceTime;
            meteoData(i).air_temperature=meteoRow.Var14 + calibrationTool.zeroDegInKelvin;
            meteoData(i).tod = hours(meteoData(i).dateTime-meteoData(1).dateTime);
            
            %rowPrec = find(isbetween(precipitation.dateTime,meteoData(i).dateTime,meteoData(i).dateTime+dt));
            %meteoData(i).precipitation = sum(precipitation(rowPrec,:).p1,'omitnan');
            meteoData(i).precipitation = nan;
            meteoData(i).rel_humidity=meteoRow.Var13;
            
            meteoData(i).air_pressure= meteoRow.Var9;
            %         % TODO Check units
            %         if i>1
            %             meteoData(i).precipitation = meteoData(i).rain_accumulation - meteoData(i-1).rain_accumulation;
            %         end
        end
    else
        error('No Exwi Meteo data !')
    end
catch ME
    warning('no meteo from local station, trying with anetz data')
    try
        meteoData = read_meteo_data_MCH(calibrationTool);
        if ~isempty(fieldnames(meteoData))
            disp('Found anetz data for this day');
        else
            error('No ANETZ data for this day !')
        end
    catch ME
        warning(ME.message)
        disp('no meteo data loaded for this day')
        meteoData = struct();
    end
end
end
