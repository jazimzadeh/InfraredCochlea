%% - Make structure w/data from all mice
clear
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets")
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";

Table       = readtable('Inventory_VGLUT3_WT.csv');
TableHeader = Table.Properties.VariableNames;

for j = 1:size(Table,1)   % Iterate through the table rows (mice); create Data structure
    Data(j).Date            = Table.Date(j);
    Data(j).Mouse           = Table.Mouse(j);
    %Data(j).Genotype_DT1_WT0= Table.Genotype_DT1_WT0(j);
    Data(j).Genotype_VGLUT_1_WT_0       = Table.Genotype_VGLUT_1_WT_0(j);
    Data(j).Strain_C57Bl6_1_CBACaJ_0    = Table.Strain_C57Bl6_1_CBACaJ_0(j);

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
save("VGLUT3_WT_Inventory.mat","Data")

%% - Load the Data structure 
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
load ("VGLUT3_WT_Inventory.mat")

% Pull out the microphonic CAP for the 90 dB click
for j = 1:size(Data,2)
    mean_all = mean(Data(j).Run_Click_90_all.*1000000,1); % Convert to uV by mult by 1e6
    mean_odd = mean(Data(j).Run_Click_90_all(2:2:end,:).*1000000,1);

    Data(j).Click90microphonic = mean_odd - mean_all;
end

% Keep only the mice that are C57Bl6
DataC57 = Data([Data.Strain_C57Bl6_1_CBACaJ_0] == 1);

% Make a structure each for the VGLUT3 and WT mice (C57Bl6 only)
DataVGLUT3  = DataC57([DataC57.Genotype_VGLUT_1_WT_0]==1);
DataWT      = DataC57([DataC57.Genotype_VGLUT_1_WT_0]==0);

% Make the time vector
% Fs = 25000;            % Sampling rate in Hz
% tstep = 1/Fs;
% t = [0:(length(DataDT(1).Run_100_200um)-49)]; % # of time points
% t = t.*tstep;           % Time vector

t  = linspace(0,0.01098,268).*1000;
tt = linspace(0,0.01098,2144).*1000;

%% - Plot all WT and all VGLUT3 mouse aCAPs, aMicrophonics, oCAPs with 400 um fiber
close all

% Settings for plot offset and y-limit 
ylimvec_aCAP = [-200 2500];
offset_aCAP  = 240;
ylimvec_oCAP = [-200 2200];
offset_oCAP  = 200;
xlimvec      = [0 10];

% Define the times within which to find aCAP amplitude
min_start   = 2.0; % Input time (in ms) at which to start looking for minimum
min_end     = 2.35;
max_start   = 2.4;
max_end     = 3.0;

% Define the times within which to find oCAP peak 1 amplitude
oCAP1_min_start   = 0.9; 
oCAP1_min_end     = 1.1;
oCAP1_max_start   = 1.1;
oCAP1_max_end     = 1.3;

% Define the times within which to find oCAP peak 2 amplitude
oCAP2_min_start   = 1.6; 
oCAP2_min_end     = 2.1;
oCAP2_max_start   = 2.1; 
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
WT_aCAP_Ampl            = zeros(size(DataWT,2),1);
VGLUT3_aCAP_Ampl        = zeros(size(DataVGLUT3,2),1);
WT_oCAP1_Ampl           = zeros(size(DataWT,2),1);
VGLUT3_oCAP1_Ampl       = zeros(size(DataVGLUT3,2),1);
WT_oCAP2_Ampl           = zeros(size(DataWT,2),1);
VGLUT3_oCAP2_Ampl       = zeros(size(DataVGLUT3,2),1);

subplot_tight(1,5,1)
offset = 0;
for p = 1:size(DataWT,2)
    signal = 0.001.*DataWT(p).Run_Click_Avg(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on    
    text(4.5,offset+8, strcat(DataWT(p).Date, ' - Ms ', num2str(DataWT(p).Mouse)),'interpreter','none')   
    offset = offset + offset_aCAP/10;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('WT - ABR 90 dB Click')

subplot_tight(1,5,2)
offset = 0;
for p = 1:size(DataWT,2)
    signal = 0.001.*DataWT(p).Run_Click_Avg(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on    
    %text(4.5,offset+80, strcat(DataWT(p).Date, ' - Ms ', num2str(DataWT(p).Mouse)),'interpreter','none')
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    WT_aCAP_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_aCAP;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP 90 dB Click')

subplot_tight(1,5,3)
offset = 0;
for p = 1:size(DataWT,2)
    plot(tt, DataWT(p).Click90microphonic(1,:) + offset, 'k','linewidth',2)
    hold on    
    %text(5,offset+80, strcat(DataWT(p).Date, ' - Ms ', num2str(DataWT(p).Mouse)),'interpreter','none')
    offset = offset + offset_aCAP;
    grid on
    %ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP Micorphonic - Click ')

subplot_tight(1,5,4)
offset = 0;
for p = 1:size(DataWT,2)
    signal = 0.001.*DataWT(p).Run_100_400um(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on     
    offset = offset + offset_aCAP/10;
    grid on   
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('ABR 100% 400 um fiber')

subplot_tight(1,5,5)
offset = 0;
for p = 1:size(DataWT,2)
    signal = 0.001.*DataWT(p).Run_100_400um(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on 
    % Find peak 1
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
    WT_oCAP1_Ampl(p) = Ampl;
    % Markers for peak 1
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Find peak 2
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
    WT_oCAP2_Ampl(p) = Ampl;
    % Markers for peak 2
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Plot 40% laser power:
    %plot(t, 0.001.*DataWT(p).Run_40_400um(3,49:316) + offset,'r','linewidth',2)
    offset = offset + offset_oCAP;
    grid on   
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP 100% 400 um fiber')

set(gcf,'position',[0 0 1500 900])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_VGLUT3_WT_AllTraces_WT.pdf'), 'ContentType', 'vector');

% - - - - - - - - - - -%
% - VGLUT3 - %
figure

subplot_tight(1,5,1)
offset = 0;
for p = 1:size(DataVGLUT3,2)
    signal = 0.001.*DataVGLUT3(p).Run_Click_Avg(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on
    text(4.5,offset+5, strcat(DataVGLUT3(p).Date, ' - Ms ', num2str(DataVGLUT3(p).Mouse)),'interpreter','none')
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    if isempty(Ampl)
        VGLUT3_aCAP_Ampl(p) = NaN;
    else
        VGLUT3_aCAP_Ampl(p) = Ampl;
    end
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_aCAP/10;
    grid on
    xlabel('ms')
    ylim(ylimvec_aCAP/10)
    xlim(xlimvec)
end
title('VGLUT3 - ABR 90 dB Click')

subplot_tight(1,5,2)
offset = 0;
for p = 1:size(DataVGLUT3,2)
    signal = 0.001.*DataVGLUT3(p).Run_Click_Avg(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on
    %text(4.5,offset+20, strcat(DataVGLUT3(p).Date, ' - Ms ', num2str(DataVGLUT3(p).Mouse)),'interpreter','none')
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    if isempty(Ampl)
        VGLUT3_aCAP_Ampl(p) = NaN;
    else
        VGLUT3_aCAP_Ampl(p) = Ampl;
    end
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_aCAP;
    grid on
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP 90 dB Click')

subplot_tight(1,5,3)
offset = 0;
for p = 1:size(DataVGLUT3,2)
    plot(tt, DataVGLUT3(p).Click90microphonic(1,:) + offset, 'k','linewidth',2)
    hold on        
    offset = offset + offset_aCAP;
    grid on
    %ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP Microphonic')

subplot_tight(1,5,4)
offset = 0;
for p = 1:size(DataVGLUT3,2)
    signal = 0.001.*DataVGLUT3(p).Run_100_400um(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on     
    offset = offset + offset_aCAP;
    grid on   
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('ABR - 400 um fiber, 100%')

subplot_tight(1,5,5)
offset = 0;
for p = 1:size(DataVGLUT3,2)
    signal = 0.001.*DataVGLUT3(p).Run_100_400um(3,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on 
    % Find peak 1
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
    if isempty(Ampl)
        VGLUT3_oCAP1_Ampl(p) = 0;
    else
        VGLUT3_oCAP1_Ampl(p) = Ampl;
    end
    % Markers for peak 1
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Find peak 2
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
    if isempty(Ampl)
        VGLUT3_oCAP2_Ampl(p) = 0;
    else
        VGLUT3_oCAP2_Ampl(p) = Ampl;
    end
    % Markers for peak 2
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Plot 40% laser power:
    %plot(t, 0.001.*DataWT(p).Run_40_400um(3,49:316) + offset,'r','linewidth',2)
    offset = offset + offset_aCAP;
    grid on   
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('CAP - 400 um fiber, 100%')

set(gcf,'position',[0 0 1500 900])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_VGLUT3_WT_AllTraces_VGLUT3.pdf'), 'ContentType', 'vector');


%% - Plot all WT and all VGLUT3 mouse aCAPs, aMicrophonics, oCAPs with 200 um fiber
clf
ylimvec_oCAP = [-200 2200];
offset_oCAP  = 200;
ylimvec_aCAP = [-200 2500];
xlimvec      = [0 10];
offset_aCAP  = 240;

subplot_tight(1,6,1)
offset = 0;
offset_step = 300;
for p = 1:size(DataWT,2)
    plot(t, 0.001.*DataWT(p).Run_Click_Avg(3,49:316) + offset, 'k','linewidth',2)
    hold on    
    text(4.5,offset+80, strcat(DataWT(p).Date, ' - Ms ', num2str(DataWT(p).Mouse)),'interpreter','none')
    offset = offset + offset_aCAP;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('WT, 90 dB Click')

subplot_tight(1,6,2)
offset = 0;
offset_step = 300;
for p = 1:size(DataWT,2)
    plot(tt, DataWT(p).Click90microphonic(1,:) + offset, 'k','linewidth',2)
    hold on    
    %text(5,offset+80, strcat(DataWT(p).Date, ' - Ms ', num2str(DataWT(p).Mouse)),'interpreter','none')
    offset = offset + offset_aCAP;
    grid on
    %ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('90 dB Click - Microphonic')

subplot_tight(1,6,3)
offset = 0;
offset_step = 120;
for p = 1:size(DataWT,2)
    if strcmp(DataWT(p).Run_100_200um, 'missing') == 1
        text(1,offset,'Missing')
    else
        plot(t, 0.001.*DataWT(p).Run_100_200um(3,49:316) + offset, 'k','linewidth',2)
    end
    hold on
    if strcmp(DataWT(p).Run_40_200um, 'missing') == 1
        text(1,offset,'Missing')
    else
        plot(t, 0.001.*DataWT(p).Run_40_200um(3,49:316) + offset, 'r','linewidth',2)
    end
    offset = offset + offset_oCAP;
    grid on   
    xlabel('ms')
    ylim(ylimvec_oCAP)
    xlim(xlimvec)
end
title('200 um fiber, 100%, 40%')

subplot_tight(1,6,4)
offset = 0;
offset_step = 270;
for p = 1:size(DataVGLUT3,2)
    plot(t, 0.001.*DataVGLUT3(p).Run_Click_Avg(3,49:316) + offset, 'k','linewidth',2)
    hold on
    %plot(t.*1000, 0.001.*DataWT(p).Run_100_200um(3,49:end) + offset, 'r','linewidth',2)
    %plot(t, 0.001.*DataWT(p).Run_40_200um(3,49:316) + offset, 'r','linewidth',2)
    text(4.5,offset+20, strcat(DataVGLUT3(p).Date, ' - Ms ', num2str(DataVGLUT3(p).Mouse)),'interpreter','none')
    offset = offset + offset_aCAP;
    grid on
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('VGLUT3, 90 dB Click')

subplot_tight(1,6,5)
offset = 0;
offset_step = 270;
for p = 1:size(DataVGLUT3,2)
    plot(tt, DataVGLUT3(p).Click90microphonic(1,:) + offset, 'k','linewidth',2)
    hold on        
    offset = offset + offset_aCAP;
    grid on
    %ylabel('uV')
    xlabel('ms')
    ylim(ylimvec_aCAP)
    xlim(xlimvec)
end
title('90 dB Click - Microphonic')

subplot_tight(1,6,6)
offset = 0;
offset_step = 100;
for p = 1:size(DataVGLUT3,2)
    if strcmp(DataVGLUT3(p).Run_100_200um, 'missing') == 1
        text(1,offset,'Missing')
    else
        plot(t, 0.001.*DataVGLUT3(p).Run_100_200um(3,49:316) + offset, 'k','linewidth',2)
    end
    hold on
    if strcmp(DataVGLUT3(p).Run_40_200um, 'missing') == 1
        text(1,offset,'Missing')
    else
        plot(t, 0.001.*DataVGLUT3(p).Run_40_200um(3,49:316) + offset, 'r','linewidth',2)  
    end
    offset = offset + offset_oCAP;
    grid on
    xlabel('ms')
    ylim(ylimvec_oCAP)
    xlim(xlimvec)
end
title('200 um fiber, 100%, 40%')

set(gcf,'position',[0 0 1500 900])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_VGLUT3_WT_oCAPs_200um.pdf'), 'ContentType', 'vector');

%% - Calculate the microphonic energy

% WT_microphonic_E        = zeros(size(DataWT,2),1);
% VGLUT3_microphonic_E    = zeros(size(DataVGLUT3,2),1);
% 
% for p = 1:size(DataWT,2)
%     signal              = DataWT(p).Click90microphonic(1,:);
%     WT_microphonic_E(p) = sum(abs(signal).^2) ; % in uV - ms    
% end
% 
% for p = 1:size(DataVGLUT3,2)
%     signal                  = DataVGLUT3(p).Click90microphonic(1,:);
%     VGLUT3_microphonic_E(p) = sum(abs(signal).^2); % in uV - ms    
% end

%% - Calculate the microphonic RMS
close all
WT_microphonic_E        = zeros(size(DataWT,2),1);
VGLUT3_microphonic_E    = zeros(size(DataVGLUT3,2),1);

for p = 1:size(DataWT,2)
    signal              = DataWT(p).Click90microphonic(1,:);
    WT_microphonic_E(p) = sqrt(sum(signal.^2));  
end

for p = 1:size(DataVGLUT3,2)
    signal                  = DataVGLUT3(p).Click90microphonic(1,:);
    VGLUT3_microphonic_E(p) = sqrt(sum(signal.^2));   
end


%% - Generate and save tables with all the summary data

SummaryWT = table(WT_aCAP_Ampl, WT_oCAP1_Ampl, WT_oCAP2_Ampl, WT_microphonic_E);
for f = 1:size(SummaryWT,1)
    SummaryWT.Date(f)     = DataWT(f).Date;
    SummaryWT.Mouse(f)    = DataWT(f).Mouse;
end

SummaryVGLUT3 = table(VGLUT3_aCAP_Ampl, VGLUT3_oCAP1_Ampl, VGLUT3_oCAP2_Ampl, VGLUT3_microphonic_E);
for f = 1:size(SummaryVGLUT3,1)
    SummaryVGLUT3.Date(f)     = DataVGLUT3(f).Date;
    SummaryVGLUT3.Mouse(f)    = DataVGLUT3(f).Mouse;
end

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
writetable(SummaryWT,'VGLUT3_WT_Summary_WT.csv')
writetable(SummaryVGLUT3,'VGLUT3_WT_Summary_VGLUT3.csv')

%% - Plot the summary data

% Test for normality
kstest(WT_aCAP_Ampl)
% All of the distributions are not normal, so we can't use parametric tests

% Levene Test for Homoscedasticity (assumption that variance is the same)
%Levenetest()
A = [[WT_aCAP_Ampl, ones(length(WT_aCAP_Ampl),1)];[VGLUT3_aCAP_Ampl, 2.*ones(length(VGLUT3_aCAP_Ampl),1)]];
A(10,:) = [];
Levenetest(A)
% This shows that the variances between the two groups are not the same

Levenetest([[WT_microphonic_E, ones(length(WT_microphonic_E),1)];[VGLUT3_microphonic_E, 2.*ones(length(VGLUT3_microphonic_E),1)]])

% Kruskal-Wallis Test - No assumptions about normality or variance - aCAP
Group1 = WT_aCAP_Ampl;
Group2 = VGLUT3_aCAP_Ampl;

A = [[Group1, ones(length(Group1),1)];[Group2, 2.*ones(length(Group2),1)]];
A(10,:) = [];
p = kruskalwallis(A(:,1),A(:,2));

% Kruskal-Wallis Test - Microphonic
Group1 = WT_microphonic_E;
Group2 = VGLUT3_microphonic_E;

A = [[Group1, ones(length(Group1),1)];[Group2, 2.*ones(length(Group2),1)]];
p = kruskalwallis(A(:,1),A(:,2));

% Kruskal-Wallis Test - oCAP1
Group1 = WT_oCAP1_Ampl;
Group2 = VGLUT3_oCAP1_Ampl;

A = [[Group1, ones(length(Group1),1)];[Group2, 2.*ones(length(Group2),1)]];
p = kruskalwallis(A(:,1),A(:,2));

% Kruskal-Wallis Test - oCAP2
Group1 = WT_oCAP2_Ampl;
Group2 = VGLUT3_oCAP2_Ampl;

A = [[Group1, ones(length(Group1),1)];[Group2, 2.*ones(length(Group2),1)]];
p = kruskalwallis(A(:,1),A(:,2));


close all
figure
subplot(2,2,1)
plot(0.4.*(rand(1,length(WT_aCAP_Ampl))),WT_aCAP_Ampl,'.k','markersize',30)
hold on
plot(0.4+0.4.*(rand(1,length(VGLUT3_aCAP_Ampl))),VGLUT3_aCAP_Ampl,'.r','markersize',30)
%[H, p] = unpaired_ttest(WT_aCAP_Ampl, VGLUT3_aCAP_Ampl, 1, 400, [0.2 0.6]);
%[p] = KruskalWallisPlot(WT_aCAP_Ampl, DT_aCAP_Ampl, 1, 400, [0.2 0.6]);
xlim([0 0.8])
ylim([0 450])
xticks([0.2 0.6])
set(gca,'fontsize',16,'xticklabel',{'WT', 'VGLUT3'})
ylabel('Peak-peak (uV)')
title('90 dB Click CAP')

subplot(2,2,2)
plot(0.7+0.4.*(rand(1,length(WT_microphonic_E))),WT_microphonic_E,'.k','markersize',40)
hold on
plot(2.7+0.5.*(rand(1,length(VGLUT3_microphonic_E))),VGLUT3_microphonic_E,'.r','markersize',40)
%[H, p] = unpaired_ttest(WT_microphonic_E, VGLUT3_microphonic_E, 1, 2000, [1 3]);
[p] = KruskalWallisPlot(WT_microphonic_E, VGLUT3_microphonic_E, 1, 2000, [1 3]);
xlim([0 4])
ylim([0 2200])
xticks([1 3])
set(gca,'fontsize',16,'xticklabel',{'WT', 'VGLUT3'})
ylabel('RMS (uV)')
title('90 dB Click CAP Microphonic')

subplot(2,2,3)
plot(0.7+0.4.*(rand(1,length(WT_oCAP1_Ampl))),WT_oCAP1_Ampl,'.k','markersize',40)
hold on
plot(2.7+0.5.*(rand(1,length(VGLUT3_oCAP1_Ampl))),VGLUT3_oCAP1_Ampl,'.r','markersize',40)
%[H, p] = unpaired_ttest(WT_oCAP1_Ampl, VGLUT3_oCAP1_Ampl, 1, 400, [1 3]);
[p] = KruskalWallisPlot(WT_oCAP1_Ampl, VGLUT3_oCAP1_Ampl, 1, 400, [1 3]);
xlim([0 4])
ylim([0 450])
xticks([1 3])
set(gca,'fontsize',16,'xticklabel',{'WT', 'VGLUT3'})
ylabel('Peak-peak (uV)')
title('oCAP1 400 um')

subplot(2,2,4)
plot(0.7+0.4.*(rand(1,length(WT_oCAP2_Ampl))),WT_oCAP2_Ampl,'.k','markersize',40)
hold on
plot(2.7+0.5.*(rand(1,length(VGLUT3_oCAP2_Ampl))),VGLUT3_oCAP2_Ampl,'.r','markersize',40)
%[H, p] = unpaired_ttest(WT_oCAP2_Ampl, VGLUT3_oCAP2_Ampl, 1, 400, [1 3]);
[p] = KruskalWallisPlot(WT_oCAP2_Ampl, VGLUT3_oCAP2_Ampl, 1, 400, [1 3]);
xlim([0 4])
ylim([0 450])
xticks([1 3])
set(gca,'fontsize',16,'xticklabel',{'WT', 'VGLUT3'})
ylabel('Peak-peak (uV)')
title('oCAP2 400 um')


set(gcf,'position',[0 100 600 500])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_VGLUT3_WT_400um_SummaryData.pdf'), 'ContentType', 'vector');

%% - Microphonic RMS vs oCAP peak 1
clf
plot(WT_microphonic_E, WT_oCAP1_Ampl,'.k','markersize',20)
hold on
plot(VGLUT3_microphonic_E, VGLUT3_oCAP1_Ampl,'.r','markersize',20)

text(10,280,'WT','fontsize',20)
text(10,250,'VGLUT3','color','r','fontsize',20)
xlabel('Microphonic RMS')
ylabel('Optical Peak 1 Amplitude')
title('Microphonic RMS vs oCAP peak 1 amplitude')

set(gca,'fontsize',16)
set(gcf,'position',[1100 600 400 300])

all_microphonic_E   = [WT_microphonic_E; VGLUT3_microphonic_E];
all_oCAP1           = [WT_oCAP1_Ampl; VGLUT3_oCAP1_Ampl];
[R, Pee] = corrcoef(all_microphonic_E, all_oCAP1);
text(1300, 50,strcat('R = ',num2str(round(R(1,2),2))),'fontsize',16)

% Save the summary numbers
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
save("VGLUT3_WT_microphonic_E_WT.mat"    ,"WT_microphonic_E")
save("VGLUT3_WT_microphonic_E_VGLUT3.mat","VGLUT3_microphonic_E")
save("VGLUT3_WT_oCAP1_Ampl_WT.mat"       ,"WT_oCAP1_Ampl")
save("VGLUT3_WT_oCAP1_Ampl_VGLUT3.mat"   ,"VGLUT3_oCAP1_Ampl")



cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_VGLUT3_WT_400um_MicrophonicVs_oCAP1.pdf'), 'ContentType', 'vector');

%% - Overlay VGLUT3 and WT laser traces - for Figure 2

clf
clear
exptdate = '2023_3_29';
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
path = strcat(folder,exptdate,'/');

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/WorkspaceVariables")
WT = importdata('2023_3_28_WT_oCAP.mat');

cd(path)
Runs_400        = [3   4   5   6   7   8  ];
PowerPercent400 = [100 80  60  40  20  17 ];
PulseDuration400= [100 100 100 100 100 100];

Target_Distance = 300;
%Pulse_Length    = 100;
Fiber_Diameter  = 400;

A = readmatrix(strcat('Run', num2str(Runs_400(1)), '.csv'));
Fs = 25000;            % Sampling rate in Hz
tstep = 1/Fs;
t = [0:(length(A)-49)]; % # of time points
t = t.*tstep;           % Time vector

cd(path)
clear u
% Make 3-D array containing all the CSV data
for i = 1:length(Runs_400)
    u(i,:,:) = readmatrix(strcat('Run', num2str(Runs_400(i)), '.csv'));
end

close all
subplot(2,1,1)
clear data EnergyDensity400 SP_400 AP_400
offset = 0;
j=1;
    % Record the power and radiant energy delivered
    [Watts400(j), ~, ~, EnergyDensity400(j)] = CapellaEnergy(Fiber_Diameter, PowerPercent400(j), Target_Distance, PulseDuration400(j));
    data = squeeze(u(j,3,49:end));

    plot(t.*1000, 0.001.*(1.0.*(data)-offset),'linewidth',2,'color',[1 0 0]) % VGLUT3
    hold on
    plot(t.*1000, 0.001.*(WT(:,j)-offset),'linewidth',2,'color',[0 0 0]) % WT

legend(strcat(num2str(EnergyDensity400(1)), 'mJ/cm^2'),'location','southeast')
set(gca,'fontsize',16)
xlim([0.5 3])
ylim([-100 55])
title('oCAPs with 400 um Fiber','Black: WT (3/28 Ms1), Red: VGLUT3 (3/29 Ms1)')

ylabel('uV')
xlabel('ms')

subplot(2,1,2)
    plot(t.*1000, 0.002.*(1.0.*(data)-offset),'linewidth',2,'color',[1 0 0]) % VGLUT3
    hold on
    plot(t.*1000, 0.001.*(WT(:,j)-offset),'linewidth',2,'color',[0 0 0]) % WT
    xlim([0.5 3])
    ylim([-100 55])
    set(gca,'fontsize',16)
    set(gcf,'Position',[900 100 600 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figure2_VGLUT3_WT_Overlaid.pdf'), 'ContentType', 'vector');


