%% - Make structure w/data from all mice
clear
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets")
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";

Table       = readtable('Inventory_TMC1_Het.csv');
TableHeader = Table.Properties.VariableNames;

for j = 1:size(Table,1)   % Iterate through the table rows (mice); create Data structure
    Data(j).Date            = Table.Date(j);
    Data(j).Mouse           = Table.Mouse(j);
    Data(j).Genotype_TMC1_Het0= Table.Genotype_TMC1_Het0(j);
    Data(j).Strain          = Table.Strain(j);

    DataPath                = cd(strcat(folder,Table.Date(j)));    
    if isnan(Table.Run_100_200um(j)) 
        Data(j).Run_100_200um    = 'missing';
    elseif ~isnan(Table.Run_100_200um(j)) 
        Data(j).Run_100_200um   = readmatrix(strcat('Run',num2str(Table.Run_100_200um(j)),'.csv'));
    end

    if isnan(Table.Run_40_200um(j)) 
        Data(j).Run_40_200um    = 'missing';
    elseif ~isnan(Table.Run_40_200um(j)) 
        Data(j).Run_40_200um    = readmatrix(strcat('Run',num2str(Table.Run_40_200um(j)),'.csv'));
    end

    if isnan(Table.Run_100_400um(j)) 
        Data(j).Run_100_400um    = 'missing';
    elseif ~isnan(Table.Run_100_400um(j)) 
        Data(j).Run_100_400um   = readmatrix(strcat('Run',num2str(Table.Run_100_400um(j)),'.csv'));
    end

    if isnan(Table.Run_40_400um(j)) 
        Data(j).Run_40_400um    = 'missing';
    elseif ~isnan(Table.Run_40_400um(j)) 
        Data(j).Run_40_400um    = readmatrix(strcat('Run',num2str(Table.Run_40_400um(j)),'.csv'));
    end
    
    Data(j).Run_Click_Avg   = readmatrix(strcat('Run',num2str(Table.Run_Click(j)),'.csv'));

    if isnan(Table.Run_eABR(j)) 
        Data(j).Run_eABR    = 'missing';
    elseif ~isnan(Table.Run_eABR(j)) 
        Data(j).Run_eABR    = readmatrix(strcat('Run',num2str(Table.Run_eABR(j)),'.csv'));
    end

    if isnan(Table.Current_eABR(j)) 
        Data(j).Current_eABR    = 'missing';
    elseif ~isnan(Table.Current_eABR(j)) 
        Data(j).Current_eABR    = Table.Current_eABR(j);
    end
       
    Data(j).Run_Click_90_all= readmatrix(strcat('Run',num2str(Table.Run_Click(j)),'-0-1-2-1.csv')); %The 90 dB CAP
end

% Save the Data structure to a mat file
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
save("TMC1_Het_Inventory.mat","Data")

%% - Load the Data structure 
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
load ("TMC1_Het_Inventory.mat")

% Pull out the microphonic CAP for the 90 dB click
for j = 1:size(Data,2)
    mean_all = mean(Data(j).Run_Click_90_all.*1000000,1); % Convert to uV by mult by 1e6
    mean_odd = mean(Data(j).Run_Click_90_all(2:2:end,:).*1000000,1);

    Data(j).Click90microphonic = mean_odd - mean_all;
end

% Make a structure each for the DT and WT mice
DataTMC1= Data([Data.Genotype_TMC1_Het0]==1);
DataHet = Data([Data.Genotype_TMC1_Het0]==0);

t  = linspace(0,0.01098,268).*1000;
tt = linspace(0,0.01098,2144).*1000;

%% - Plot all WT and all DT mouse aCAPs, aMicrophonics, oCAPs
clf
% Settings for plot offset and y-limit 
ylimvec_aCAP = [-200 800];
offset_aCAP  = 240;
ylimvec_oCAP = [-200 800];
offset_oCAP  = 240;
xlimvec      = [0 10];

% Define the times within which to find aCAP amplitude
min_start   = 1.8; % Input time (in ms) at which to start looking for minimum
min_end     = 2.6;
max_start   = 2.35;
max_end     = 3.2;

% Define the times within which to find oCAP peak 1 amplitude
oCAP1_min_start   = 0.9; 
oCAP1_min_end     = 1.1;
oCAP1_max_start   = 1.1;
oCAP1_max_end     = 1.4;

% Define the times within which to find oCAP peak 2 amplitude
oCAP2_min_start   = 1.6; 
oCAP2_min_end     = 2.1;
oCAP2_max_start   = 2.0; 
oCAP2_max_end     = 2.5;

% Convert the times in ms to indices of the vector t
[~, min_start_idx]  = findNearest(t, min_start);
[~, min_end_idx]    = findNearest(t, min_end);
[~, max_start_idx]  = findNearest(t, max_start);
[~, max_end_idx]    = findNearest(t, max_end);

% Convert the times in ms to indices of the vector t
[~, oCAP1_min_start_idx]  = findNearest(t, oCAP1_min_start);
[~, oCAP1_min_end_idx]    = findNearest(t, oCAP1_min_end);
[~, oCAP1_max_start_idx]  = findNearest(t, oCAP1_max_start);
[~, oCAP1_max_end_idx]    = findNearest(t, oCAP1_max_end);

% Convert the times in ms to indices of the vector t
[~, oCAP2_min_start_idx]  = findNearest(t, oCAP2_min_start);
[~, oCAP2_min_end_idx]    = findNearest(t, oCAP2_min_end);
[~, oCAP2_max_start_idx]  = findNearest(t, oCAP2_max_start);
[~, oCAP2_max_end_idx]    = findNearest(t, oCAP2_max_end);

% Initialize the summary variables
Het_aCAP_Ampl           = zeros(size(DataHet,2),1);
TMC1_aCAP_Ampl          = zeros(size(DataTMC1,2),1);
Het_oCAP1_Ampl          = zeros(size(DataHet,2),1);
TMC1_oCAP1_Ampl         = zeros(size(DataTMC1,2),1);
Het_oCAP2_Ampl          = zeros(size(DataHet,2),1);
TMC1_oCAP2_Ampl         = zeros(size(DataTMC1,2),1);

subplot_tight(1,5,1)
offset = 0;
for p = 1:size(DataHet,2)
    signal = 0.001.*DataHet(p).Run_Click_Avg(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on    
    text(4.5,offset+6, strcat(DataHet(p).Date, ' - Ms ', num2str(DataHet(p).Mouse)),'interpreter','none')
    offset = offset + offset_aCAP/10;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('WT (Het) ABR, 90 dB Click')

subplot_tight(1,5,2)
offset = 0;
for p = 1:size(DataHet,2)
    signal = 0.001.*DataHet(p).Run_Click_Avg(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on    
    text(4.5,offset+60, strcat(DataHet(p).Date, ' - Ms ', num2str(DataHet(p).Mouse)),'interpreter','none')
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    Het_aCAP_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_aCAP;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('WT (Het) CAP, 90 dB Click')

subplot_tight(1,5,3)
offset = 0;
for p = 1:size(DataHet,2)
    plot(tt, DataHet(p).Click90microphonic(1,:) + offset, 'k','linewidth',2)
    hold on    
    %text(5,offset+80, strcat(DataHet(p).Date, ' - Ms ', num2str(DataHet(p).Mouse)),'interpreter','none')
    offset = offset + offset_aCAP;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP Click - Microphonic')

subplot_tight(1,5,4)
offset = 0;
for p = 1:size(DataHet,2)
    signal = 0.001.*DataHet(p).Run_100_400um(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on  
    offset = offset + offset_aCAP/10;
    grid on   
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('ABR 400 um fiber, 100%')

subplot_tight(1,5,5)
offset = 0;
for p = 1:size(DataHet,2)
    signal = 0.001.*DataHet(p).Run_100_400um(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on  
     % Find peak 1
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
    Het_oCAP1_Ampl(p) = Ampl;
    % Markers for peak 1
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Find peak 2
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
    Het_oCAP2_Ampl(p) = Ampl;
    % Markers for peak 2
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Plot 40% power, 400 um
    %plot(t, 0.001.*DataHet(p).Run_40_400um(3,49:316) + offset, 'r','linewiTMC1h',2)
    offset = offset + offset_oCAP;
    grid on   
    xlabel('ms')
    ylim(ylimvec_oCAP)
    xlim(xlimvec)
end
title('CAP 400 um fiber, 100%')

set(gcf,'position',[0 0 1500 900])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_TMC1_Het_AllTraces_Het.pdf'), 'ContentType', 'vector');


% - - - - - - - - - - %
% - TMC1 - %
figure
subplot_tight(1,5,1)
offset = 0;
for p = 1:size(DataTMC1,2)
    signal = 0.001.*DataTMC1(p).Run_Click_Avg(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on    
    offset = offset + offset_aCAP/10;
    grid on
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('ABR TMC1, 90 dB Click')

subplot_tight(1,5,2)
offset = 0;
for p = 1:size(DataTMC1,2)
    signal = 0.001.*DataTMC1(p).Run_Click_Avg(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on
    %plot(t.*1000, 0.001.*DataHet(p).Run_100_200um(3,49:end) + offset, 'r','linewiTMC1h',2)
    %plot(t, 0.001.*DataHet(p).Run_40_400um(3,49:316) + offset, 'r','linewiTMC1h',2)
    text(4.5,offset+20, strcat(DataTMC1(p).Date, ' - Ms ', num2str(DataTMC1(p).Mouse)),'interpreter','none')
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    if isempty(Ampl)
        TMC1_aCAP_Ampl(p) = NaN;
    else
        TMC1_aCAP_Ampl(p) = Ampl;
    end
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_aCAP;
    grid on
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP TMC1, 90 dB Click')

subplot_tight(1,5,3)
offset = 0;
for p = 1:size(DataTMC1,2)
    plot(tt, DataTMC1(p).Click90microphonic(1,:) + offset, 'k','linewidth',2)
    hold on        
    offset = offset + offset_aCAP;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP Click - Microphonic')

subplot_tight(1,5,4)
offset = 0;
%offset_step = 100;
for p = 1:size(DataTMC1,2)
    signal = 0.001.*DataTMC1(p).Run_100_400um(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on  
    offset = offset + offset_aCAP/10;
    grid on
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('ABR 400 um fiber, 100%')

subplot_tight(1,5,5)
offset = 0;
offset_step = 100;
for p = 1:size(DataTMC1,2)
    signal = 0.001.*DataTMC1(p).Run_100_400um(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on
    % Find peak 1
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
    if isempty(Ampl)
        TMC1_oCAP1_Ampl(p) = 0;
    else
        TMC1_oCAP1_Ampl(p) = Ampl;
    end
    % Markers for peak 1
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Find peak 2
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
    if isempty(Ampl)
        TMC1_oCAP2_Ampl(p) = 0;
    else
        TMC1_oCAP2_Ampl(p) = Ampl;
    end
    % Markers for peak 2
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    %plot(t, 0.001.*DataTMC1(p).Run_40_400um(3,49:316) + offset, 'r','linewiTMC1h',2)    
    offset = offset + offset_oCAP;
    grid on
    xlabel('ms')
    ylim(ylimvec_oCAP)
    xlim(xlimvec)
end
title('CAP 400 um fiber, 100%')

set(gcf,'position',[0 0 1500 900])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_TMC1_Het_AllTraces_TMC1.pdf'), 'ContentType', 'vector');

%% - Calculate the microphonic energy

Het_microphonic_E        = zeros(size(DataHet,2),1);
TMC1_microphonic_E       = zeros(size(DataTMC1,2),1);

for p = 1:size(DataHet,2)
    signal              = DataHet(p).Click90microphonic(1,:);
    Het_microphonic_E(p) = sum(abs(signal)) .* tt(2); % in uV - ms    
end

for p = 1:size(DataTMC1,2)
    signal                  = DataTMC1(p).Click90microphonic(1,:);
    TMC1_microphonic_E(p)     = sum(abs(signal)) .* tt(2); % in uV - ms    
end

%% - Calculate the microphonic RMS

Het_microphonic_E        = zeros(size(DataHet,2),1);
TMC1_microphonic_E       = zeros(size(DataTMC1,2),1);

for p = 1:size(DataHet,2)
    signal                  = DataHet(p).Click90microphonic(1,:);
    Het_microphonic_E(p)    = sqrt(sum(signal.^2));
end

for p = 1:size(DataTMC1,2)
    signal                  = DataTMC1(p).Click90microphonic(1,:);
    TMC1_microphonic_E(p)   = sqrt(sum(signal.^2));  
end

%% - Generate and save tables with all the summary data

SummaryHet = table(Het_aCAP_Ampl, Het_oCAP1_Ampl, Het_oCAP2_Ampl, Het_microphonic_E);
for f = 1:size(SummaryHet,1)
    SummaryHet.Date(f)     = DataHet(f).Date;
    SummaryHet.Mouse(f)    = DataHet(f).Mouse;
end

SummaryTMC1 = table(TMC1_aCAP_Ampl, TMC1_oCAP1_Ampl, TMC1_oCAP2_Ampl, TMC1_microphonic_E);
for f = 1:size(SummaryTMC1,1)
    SummaryTMC1.Date(f)     = DataTMC1(f).Date;
    SummaryTMC1.Mouse(f)    = DataTMC1(f).Mouse;
end

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
writetable(SummaryHet,'TMC1_Het_Summary_Het.csv')
writetable(SummaryTMC1,'TMC1_Het_Summary_TMC1.csv')

%% - Plot the summary data: action potential sizes, microphonic E
close all
subplot(2,2,1)
plot(0.4.*(rand(1,length(Het_aCAP_Ampl))),Het_aCAP_Ampl,'.k','markersize',40)
hold on
plot(0.4+0.35.*(rand(1,length(TMC1_aCAP_Ampl))),TMC1_aCAP_Ampl,'.r','markersize',40)
xlim([0 0.8])
ylim([-20 450])
xticks([0.2 0.6])
title('90 dB Click CAP')
ylabel('Peak-peak (uV)')
set(gca,'fontsize',16,'xticklabel',{'Het', 'TMC1'})
%[H, p] = unpaired_ttest(Het_aCAP_Ampl, TMC1_aCAP_Ampl, 1, 400, [0.2 0.6]);
[p] = KruskalWallisPlot(Het_aCAP_Ampl, TMC1_aCAP_Ampl, 1, 400, [0.2 0.6]);

subplot(2,2,2)
plot(0.7+0.4.*(rand(1,length(Het_microphonic_E))),Het_microphonic_E,'.k','markersize',40)
hold on
plot(2.7+0.5.*(rand(1,length(TMC1_microphonic_E))),TMC1_microphonic_E,'.r','markersize',40)
xlim([0 4])
ylim([-50 1600])
xticks([1 3])
title('90dB Click CAP Microphonic')
ylabel('RMS (uV)')
set(gca,'fontsize',16,'xticklabel',{'Het', 'TMC1'})
%[H, p] = unpaired_ttest(Het_microphonic_E, TMC1_microphonic_E, 1, 1400, [1 3]);
[p] = KruskalWallisPlot(Het_microphonic_E, TMC1_microphonic_E, 1, 1400, [1 3]);

subplot(2,2,3)
plot(0.7+0.4.*(rand(1,length(Het_oCAP1_Ampl))),Het_oCAP1_Ampl,'.k','markersize',40)
hold on
plot(2.7+0.5.*(rand(1,length(TMC1_oCAP1_Ampl))),TMC1_oCAP1_Ampl,'.r','markersize',40)
%[H, p] = unpaired_ttest(Het_oCAP1_Ampl, TMC1_oCAP1_Ampl, 1, 350, [1 3]);
[p] = KruskalWallisPlot(Het_oCAP1_Ampl, TMC1_oCAP1_Ampl, 1, 350, [1 3]);
xlim([0 4])
xticks([1 3])
ylim([-20 450])
title('oCAP1 400um')
ylabel('Peak-peak (uV)')
set(gca,'fontsize',16,'xticklabel',{'Het', 'TMC1'})


subplot(2,2,4)
plot(0.7+0.4.*(rand(1,length(Het_oCAP2_Ampl))),Het_oCAP2_Ampl,'.k','markersize',40)
hold on
plot(2.7+0.5.*(rand(1,length(TMC1_oCAP2_Ampl))),TMC1_oCAP2_Ampl,'.r','markersize',40)
xlim([0 4])
ylim([-10 250])
xticks([1 3])
title('oCAP2 400um')
ylabel('Peak-peak (uV)')
set(gca,'fontsize',16,'xticklabel',{'Het', 'TMC1'})
%[H, p] = unpaired_ttest(Het_oCAP2_Ampl, TMC1_oCAP2_Ampl, 1, 200, [1 3]);
[p] = KruskalWallisPlot(Het_oCAP2_Ampl, TMC1_oCAP2_Ampl, 1, 200, [1 3]);

set(gcf,'position',[0 100 600 500])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_TMC1_Het_400um_SummaryData.pdf'), 'ContentType', 'vector');

%% - Microphonic RMS vs oCAP peak 1
clf
plot(Het_microphonic_E, Het_oCAP1_Ampl,'.k','markersize',40)
hold on
plot(TMC1_microphonic_E, TMC1_oCAP1_Ampl,'.r','markersize',40)
%xlim([0 1400])
%ylim([0 350])
text(20,330,'WT (TMC1 +/-)','fontsize',20)
text(20,300,'TMC1 -/-','color','r','fontsize',20)
xlabel('Microphonic RMS')
ylabel('Optical Peak 1 Amplitude')
title('Microphonic RMS vs oCAP peak 1 amplitude')
ylim([-20 350])
line([0 1400], [0 0],'color',[0.8 0.8 0.8])
set(gca,'fontsize',16)
set(gcf,'position',[1000 500 500 400])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_TMC1_Het_400um_MicrophonicVs_oCAP1.pdf'), 'ContentType', 'vector');

%% - Microphonic RMS vs oCAP peak 2
clf
plot(Het_microphonic_E, Het_oCAP2_Ampl,'.k','markersize',40)
hold on
plot(TMC1_microphonic_E, TMC1_oCAP2_Ampl,'.r','markersize',40)
%xlim([0 160])
%ylim([0 250])
text(5,240,'WT (TMC1 +/-)','fontsize',20)
text(5,220,'TMC1 -/-','color','r','fontsize',20)
xlabel('Microphonic Energy')
ylabel('Optical Peak 2 Amplitude')
title('Microphonic RMS vs oCAP peak 2 amplitude')

set(gca,'fontsize',16)
set(gcf,'position',[600 500 700 600])

%% - oCAP peak 1 vs oCAP peak 2
clf
plot(Het_oCAP1_Ampl, Het_oCAP2_Ampl,'.k','markersize',40)
hold on
plot(TMC1_oCAP1_Ampl, TMC1_oCAP2_Ampl,'.r','markersize',40)
xlim([0 350])
ylim([0 200])
text(5,190,'WT (TMC1 +/-)','fontsize',20)
text(5,175,'TMC1 -/-','color','r','fontsize',20)
xlabel('Optical Peak 1 Amplitude')
ylabel('Optical Peak 2 Amplitude')
title('Microphonic Energy vs oCAP peak 1 amplitude')

set(gca,'fontsize',16)
set(gcf,'position',[600 500 700 600])


