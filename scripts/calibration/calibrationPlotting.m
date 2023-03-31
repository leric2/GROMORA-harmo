function calibrationPlotting(calibratedSpectra,calibrationTool)
%%%%% functions to plot the calibrated spectra and other %%%%
%%%%% variables from the calibration routine             %%%%

%Define plotting variables
xlims_f_wide = [-.5;.5]; %BW to plot
xlims_f_narrow = [-.2;.2]; %BW to plot

%spectrometer running mean averaging factor
specAvF = [1,5,5];

[year,month,day] = deal(string(calibratedSpectra(1).year),string(calibratedSpectra(1).month),string(calibratedSpectra(1).day));
savefile = fullfile(calibrationTool.plot_directory, year,month,day);

settings_metadata = calibrationTool.plotSettingsMetaData;

%%%%%%%  Define some variables for handling number of spectrometers  %%%%%%

try
    mkdir(savefile)
catch
    'Directory already exists'
end

%if sepectrometerQuantity not set, assume 1
if ~isfield(calibrationTool, 'spectrometerQuantity')
    calibrationTool.spectrometerQuantity = 1;
    calibrationTool.spectroSubChan = [calibrationTool.numberOfChannels];
end

initChan = 1;
for j = 1:calibrationTool.spectrometerQuantity
    stChan(j) = initChan;
    endChan(j)= stChan(j)+calibrationTool.spectroSubChan(j)-1;
    initChan = endChan(j)+1;
end

%%%%%%%%%%%%%%%%%  Plotting info in text file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save a text file containing important plot information
%this information should be saved as fields in the settings metadata struct

fn = fieldnames(settings_metadata);
comment = '';
for k=1:numel(fn)
    temp = fn{k};
    comment = append(comment, temp,': ', settings_metadata.(fn{k}), sprintf('\n'));
end

fid = fopen(append(savefile,'plot_info.txt'),'wt');
fprintf(fid, comment);
fclose(fid);

%%%%%%%%%%%%%%%% Extract Variables to plot from Structures %%%%%%%%%%%%%%%
%define the spectra variables to plot
obsFreq = vertcat(calibratedSpectra.observationFreq);%central freq
freq = vertcat(calibratedSpectra.freq); %frequency vector
plotFreq = freq/1e9 - obsFreq;

y_obs =  vertcat(calibratedSpectra.Tb);
T_sys =  vertcat(calibratedSpectra.TSys);
T_line = vertcat(calibratedSpectra.Tb_line);

%Brightness Temperauture Differences
hotDiff =  vertcat(calibratedSpectra.THotDiff); %frequency vectorTHotDiff
refDiff =  vertcat(calibratedSpectra.Tb_ref_diff); %frequency vectorTHotDiff
coldDiff =  vertcat(calibratedSpectra.TColdDiff); %frequency vectorTHotDiff

%time variables for each spectra
timeMean = vertcat(calibratedSpectra.timeMean);
dateTime = datetime(1970,1,1) + seconds(timeMean);
timeOfDay = arrayfun(@(x) timeofday(x), dateTime);

%create date subfolder (might change later to plot type)

%%%%%%%%%%%%   Defining structure for plot parameters   %%%%%%%%%%%%%  

%define struct for freq plots 
PB_freq.xlab = 'Frequency (GHz)';
PB_freq.ylab =  'Brightness Temperature (K)';
PB_freq.logy = false;
PB_freq.ydec = false;
PB_freq.xlims = xlims_f_wide;

%clear PB

%y error plots 
PB = PB_freq;
PB.multi_x = true;
PB.multi_y = true;
PB.logy = false;


for i = 1:length(stChan)
   labels(i) = "Spectrometer" + string(i);
end
 PB.labels = labels;

%plot
for i = 1:length(timeMean)
    clear PB.savestring
    PB.savestring = fullfile(savefile, "sys_T_" + string(timeOfDay(i)) + "_UTC");
    y_vector = cell(1, length(stChan));
    f_vector = cell(1, length(stChan));
    PB.ylims =[mean(T_sys(i,:))-40, mean(T_sys(i,:))+20];
    for k = 1:length(stChan)
        f_vector{k} = plotFreq(i,stChan(k):endChan(k));
        y_vector{k} = T_sys(i,stChan(k):endChan(k));
    end
    plot_basic(f_vector, y_vector, PB);
end
    

%T line plots 
%plot


PB.multi_x = true;
PB.multi_y = true;
PB.logy = false;
PB.xlims =  xlims_f_narrow;
for i = 1:length(timeMean)
    clear PB.savestring
    PB.savestring = fullfile(savefile, "Tb_line_" + string(timeOfDay(i)) + "_UTC");
    y_vector = cell(1, length(stChan));
    f_vector = cell(1, length(stChan));
    PB.ylims =[mean(T_line(i,:))-1, mean(T_line(i,:))+1];
    for k = 1:length(stChan)
        f_vector{k} = plotFreq(i,stChan(k):endChan(k));
        y_vector{k} = T_line(i,stChan(k):endChan(k));
    end
    plot_basic(f_vector, y_vector, PB);
end


%residual plots 
PB.multi_x = true;
PB.multi_y = true;
PB.logy = false;
PB.xlims =  xlims_f_narrow;
for i = 1:length(timeMean)
    clear PB.savestring
    PB.savestring = fullfile(savefile, "Tb_bal" + string(timeOfDay(i)) + "_UTC");
    y_vector = cell(1, length(stChan));
    f_vector = cell(1, length(stChan));
    PB.ylims =[mean(y_obs(i,:))-.75, mean(y_obs(i,:))+.75];
    for k = 1:length(stChan)
        f_vector{k} = plotFreq(i,stChan(k):endChan(k));
        y_vector{k} = y_obs(i,stChan(k):endChan(k));
    end
    plot_basic(f_vector, y_vector, PB);
end

%Difference Plots
PB.xlims = xlims_f_wide;
for i = 1:length(timeMean)
    clear PB.savestring
    PB.savestring = fullfile(savefile, "Tb_ref_diff" + string(timeOfDay(i)) + "_UTC");
    y_vector = cell(1, length(stChan));
    f_vector = cell(1, length(stChan));
    PB.ylims =[mean(refDiff(i,:))-.2, mean(refDiff(i,:))+.2];
    for k = 1:length(stChan)
        f_vector{k} = plotFreq(i,stChan(k):endChan(k));
        y_vector{k} = movmean(refDiff(i,stChan(k):endChan(k)),specAvF(k)* 16);
    end
    plot_basic(f_vector, y_vector, PB);
end

%Difference Plots
PB.xlims = xlims_f_wide;
for i = 1:length(timeMean)
    clear PB.savestring
    PB.savestring = fullfile(savefile, "Tb_cold_diff" + string(timeOfDay(i)) + "_UTC");
    y_vector = cell(1, length(stChan));
    f_vector = cell(1, length(stChan));
    PB.ylims =[mean(coldDiff(i,:))-.2, mean(coldDiff(i,:))+.2];
    for k = 1:length(stChan)
        f_vector{k} = plotFreq(i,stChan(k):endChan(k));
        y_vector{k} = movmean(coldDiff(i,stChan(k):endChan(k)), specAvF(k)*64);
    end
    plot_basic(f_vector, y_vector, PB);
end

PB.xlims = xlims_f_wide;
for i = 1:length(timeMean)
    clear PB.savestring
    PB.savestring = fullfile(savefile, "Tb_hot_diff" + string(timeOfDay(i)) + "_UTC");
    y_vector = cell(1, length(stChan));
    f_vector = cell(1, length(stChan));
    PB.ylims =[mean(hotDiff(i,:))-.5, mean(hotDiff(i,:))+.5];
    for k = 1:length(stChan)
        f_vector{k} = plotFreq(i,stChan(k):endChan(k));
        y_vector{k} =  movmean(hotDiff(i,stChan(k):endChan(k)),specAvF(k)*64);
    end
    plot_basic(f_vector, y_vector, PB);
end



function plot_basic(x,y,PB)
    %expected fields
    pb_fields = ["xlab" "ylab" "savestring" "logy" "ydec" "labels" "xlims" "ylims" "multi_x" "multi_y"];
    
    %set expected fields to false if not found
    for i = 1:numel(pb_fields)
      if ~ isfield(PB, pb_fields{i})
        PB.(pb_fields{i}) = false;
      end
    end
    %initialise fig
    fig = figure('visible','off');
    hold on
    grid on

    %iterate through lines to plot
    for i = 1:length(PB.labels)
       if PB.multi_y == 1
         y_temp = y{i};
       else
         y_temp = y;
       end
       if PB.multi_x == 1       
         x_temp = x{i};
       else
         x_temp = x;
       end
       plot(x_temp,y_temp,'LineWidth',3);
    end

    
    %set axes variables
    if PB.logy == 1; set(gca, 'YScale', 'log'); end
    if PB.ydec == 1; set( gca, 'YDir', 'reverse' ); end
    if ~all(PB.xlims == 0); xlim(PB.xlims); end
    if ~all(PB.ylims == 0); ylim(PB.ylims); end
    %set label details
    xlabel(PB.xlab);
    ylabel(PB.ylab);
    legend(PB.labels);
    set(gca,'FontSize',16)
    set(gca,'LineWidth',.5)

    %want to add meta-data to each fig, not sure about a good way to do
    %this
    saveas(fig,append(PB.savestring,'.png'), 'png' )
    close all
 end
 
 end
