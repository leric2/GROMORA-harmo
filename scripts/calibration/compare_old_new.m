figure()

set(gcf,'color','w');
set(gcf, 'Position',  [0, 0, 800, 600])
%plot(drift.dailyMedianHotSpectra,'b','LineWidth',1)
plot(level1b_filter.integration(3).if(level1b_filter.integration(3).potentialBadChannels==0), level1b_filter.integration(3).Tb(level1b_filter.integration(3).potentialBadChannels==0),'b','LineWidth',0.1)
hold on
plot(level1b_nofilter.integration(3).if(level1b_nofilter.integration(3).potentialBadChannels==0), level1b_nofilter.integration(3).Tb(level1b_filter.integration(3).potentialBadChannels==0),'r','LineWidth',0.1)

grid on
legend('flags filtering', 'no flag filtering', 'FontSize',12)
xlabel('IF [MHz]', 'FontSize',14)
xlim([0,1000])
ylim([60,85])
ylabel('T_B [K]', 'FontSize',14)
title('MOPI5, AC240 integrated spectra from 16-Feb-2019 12:00-18:00', 'FontSize',14,'HorizontalAlignment', 'center')

print([calibrationTool.level1Folder 'calSpectra_outleir_TN' calibrationTool.dateStr],'-dpdf')


%%
% plotting outlier examples
figure()
ind_outlier_hot=11348;

set(gcf,'color','w');
set(gcf, 'Position',  [0, 0, 800, 600])
%plot(drift.dailyMedianHotSpectra,'b','LineWidth',1)
plot(rawSpectra(ind_outlier_hot,:),'b','LineWidth',0.1)
hold on
plot(drift.dailyMedianHotSpectra+3*drift.dailyStdHotSpectra,'r','LineWidth',1)

plot(drift.dailyMedianHotSpectra-3*drift.dailyStdHotSpectra,'r','LineWidth',1)

grid on
legend('outlier, 21-Feb-2019 11:23:02', 'daily median\pm3\sigma', 'FontSize',12)
xlabel('Channels [-]', 'FontSize',14)
xlim([0,16384])
ylim([0,2.2*1e4])
ylabel('FFT counts [-]', 'FontSize',14)
%title('SOMORA, hot spectrum outlier from 21-Feb-2019', 'FontSize',16)

print([calibrationTool.level1Folder 'outlier_hot_' calibrationTool.dateStr],'-dpdf')

%%
% plotting script, old and new routine
load('/home/eric/Documents/PhD/GROSOM/Level1/level1_dt_60min_2019_10_03.mat')
somora_old_filename = '/home/eric/Documents/PhD/GROSOM/Level1/m20191003.bin';
somora_old_filename_log = '/home/eric/Documents/PhD/GROSOM/Level1/r20191003.bin';

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C'
instrumentName='SOMORA';

% Type of calibration to do: standard or debug
calibrationType='standard';

dateStr='2019_10_03';

calibrationTool=import_default_calibrationTool(instrumentName,dateStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the calibrationTool for this particular day and instrument:
% calibrationTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};

% for gaining time.
%calibrationTool.level1aExist=false;

calibrationTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
calibrationTool.calendar='standard';
labviewLog = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instrument specific parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On the long term this should be all taken from import_default_calTool

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.hotSpectraNumberOfStdDev=3;
calibrationTool.coldSpectraNumberOfStdDev=3;

calibrationTool.labviewLog = labviewLog;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug mode and plot options
calibrationTool.calType=calibrationType;

calibrationTool.rawSpectraPlot = false;
calibrationTool.calibratedSpectraPlot = true;
calibrationTool.integratedSpectraPlot = true;
% Time interval for doing the calibration
calibrationTool.calibrationTime=10;

% Total integration time
calibrationTool.integrationTime=30;
calibrationTool.minNumberOfAvgSpectra = 2;

calibrationTool.filterByTransmittance = true;
calibrationTool.filterByFlags = true;

% Temperature of the cold load
calibrationTool.TCold=80;

[calibrationTool, level1b_somora] = run_integration(calibrationTool);

% 'GROMOS' // 'SOMORA' // 'mopi5' // 'MIAWARA-C'
instrumentName='GROMOS';

% Type of calibration to do: standard or debug
calibrationType='standard';

dateStr='2019_10_03';

calibrationTool=import_default_calibrationTool(instrumentName,dateStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editing the calibrationTool for this particular day and instrument:
% calibrationTool.requiredFields={'instrumentName','bytesPerValue','rawFileFolder'};

% for gaining time.
%calibrationTool.level1aExist=false;

calibrationTool.meanDatetimeUnit='days since 1970-01-01 00:00:00';
calibrationTool.calendar='standard';
labviewLog = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instrument specific parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On the long term this should be all taken from import_default_calTool

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationTool.hotSpectraNumberOfStdDev=3;
calibrationTool.coldSpectraNumberOfStdDev=3;

calibrationTool.labviewLog = labviewLog;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug mode and plot options
calibrationTool.calType=calibrationType;

calibrationTool.rawSpectraPlot = false;
calibrationTool.calibratedSpectraPlot = true;
calibrationTool.integratedSpectraPlot = true;
% Time interval for doing the calibration
calibrationTool.calibrationTime=10;

% Total integration time
calibrationTool.integrationTime=60;
calibrationTool.minNumberOfAvgSpectra = 3;

calibrationTool.filterByTransmittance = true;
calibrationTool.filterByFlags = true;

% Temperature of the cold load
calibrationTool.TCold=80;

[calibrationTool, level1b_gromos] = run_integration(calibrationTool);

fid = fopen(somora_old_filename,'r','ieee-le');
if fid >= 3
    disp(['File ',somora_old_filename,' loaded']);
    for i=1:48
        record(i).BT         = fread( fid, 16384, 'float32' ); % [K]
        record(i).sigma_BT   = fread( fid, 16384, 'float32' ); % [K]
    end
    fclose(fid);
elseif fid == -1
    disp(['Couldn''t open file ', somora_old_filename]);
    data = [];
end

center=1.4217504e+11;
somora_frequencies=center-(8191*(1e+9/16384)):(1e+9/16384):center+(8192*(1e+9/16384));

%%
i = 38;
fig=figure()
ax = subplot(1,2,1);
set(gcf, 'Position',  [0, 0, 1100, 600])
orient(fig,'landscape')
set(gcf,'color','w');
plot(somora_frequencies(105:end)/1e9,record(i).BT(105:end),'b','LineWidth',0.1)
hold on
plot(level1b_somora.integration(i).freq(level1b_somora.integration(i).potentialBadChannels==0)/1e9,level1b_somora.integration(i).Tb(level1b_somora.integration(i).potentialBadChannels==0),'r','LineWidth',0.1)

grid on
legend('old','new','Location','northeast')
xlabel('frequency [GHz]')
%xlim([141.6,142.8])
%ylim([median(record(i).BT(105:end))-10,median(record(i).BT(105:end))+25])
ylim([95,130])
ylabel('T_B [K]')
title(['SOMORA : ',datestr(level1b_somora.integration(i).dateTime)])

i=i/2;
ax2 = subplot(1,2,2);

plot(ax2, level1b_gromos.integration(i).freq(level1b_gromos.integration(i).potentialBadChannels==0)/1e9,level1b_gromos.integration(i).Tb(level1b_gromos.integration(i).potentialBadChannels==0),'r','LineWidth',0.1)
hold on
plot(ax2, level1_data(i).frequencies/1e9,level1_data(i).bt,'b','LineWidth',0.6)
grid on
%legend('new','old')
xlabel('frequency [GHz]')
%xlim([141.6,142.8])
ylim([95,130])
%ylim([median(level1_data(i).bt)-10,median(level1_data(i).bt)+25])
ylabel('T_B [K]')
title(['GROMOS : ',datestr(level1b_gromos.integration(i).dateTime)])

print([calibrationTool.level1Folder 'comparison_GROMOS_SOMORA_77K' calibrationTool.dateStr],'-dpdf')
