%% - Fig Part B, option 1
% % --- ENTER REQUIRED INPUTS --- %%%%
% - Set the address of the mouse folder: 
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
exptdate = '011723';
mouse = 'ms2';

% - Enumerate the file #s that I wish to plot. Use the Excel file
% with experiment notes to decide which experiments to plot.
filesToPlot = 13;

% - Select plotting parameters
%t_start  = 29;      % Time in ms at which to start plotting the laser OCT and CAP
%t_end    = 34;      % Time in ms at which to stop plotting the laser OCT and CAP
%laser    = 30;      % Time at which laser was fired
laser    = 30;
t_start  = laser-1;
t_end    = laser + 4;
saveyn   = 0;       % 1 = save the plot in \Results\Oghalai\
genotype = 'WT';    % Enter mouse genotype (to properly label the plot and the file)

%%%% ------------------------------ %%%%

% - Enter the mouse folder
path = fullfile(folder,exptdate,mouse);
cd(path);

% - List all experiment files in the mouse folder
files = dir(path); files = {files.name}; files = files(4:end);

% - Loop over file to plot, importing their data and then plotting it
for f = 1 : length(filesToPlot)
    filename = char(files(filesToPlot(f)));

    [t_OCT, OCT, CAPs, t_CAP, Ampl] = SoundLaser_OCTCAP(path, filename, 'start','end',0);

%    Plot_OCTCAP_local(path, filename, mouse, OCT, t_OCT, CAPs, t_CAP, Ampl, t_start, t_end, laser, genotype, saveyn)
end

%% - Make plot
subplot(3,1,1)
plot(t_OCT,OCT{6,2},'k','linewidth',1.5)
text(19.6, 95, '50 dB')
xlim([19.5 22])
%ylim([50 100])
title(filename,'interpreter','none')
ylabel('nm')

subplot(3,1,2)
plot(t_OCT,OCT{1,2},'k','linewidth',1.5)
xlim([29.5 32])
ylim([-70 -20])

subplot(3,1,3)
plot(t_OCT,OCT{6,2},'k','linewidth',1.5)
text(29.6, 95, '50 dB')
xlim([29.5 32])
ylim([50 100])
xlabel('ms')

%% - Fig Part B, 2nd option (better) 
% - %%%% --- ENTER REQUIRED INPUTS --- %%%%
% - Set the address of the mouse folder: 
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
exptdate = '011523';
mouse = 'ms1';

% - Enumerate the file #s that I wish to plot. Use the Excel file
% with experiment notes to decide which experiments to plot.
filesToPlot = 17;%11;

% - Select plotting parameters
%t_start  = 29;      % Time in ms at which to start plotting the laser OCT and CAP
%t_end    = 34;      % Time in ms at which to stop plotting the laser OCT and CAP
%laser    = 30;      % Time at which laser was fired
laser    = 30;
t_start  = laser-1;
t_end    = laser + 4;
saveyn   = 0;       % 1 = save the plot in \Results\Oghalai\
genotype = 'WT';    % Enter mouse genotype (to properly label the plot and the file)

%%%% ------------------------------ %%%%

% - Enter the mouse folder
path = fullfile(folder,exptdate,mouse);
cd(path);

% - List all experiment files in the mouse folder
files = dir(path); files = {files.name}; files = files(4:end);

% - Loop over file to plot, importing their data and then plotting it
for f = 1 : length(filesToPlot)
    filename = char(files(filesToPlot(f)));

    [t_OCT, OCT, CAPs, t_CAP, Ampl] = SoundLaser_OCTCAP(path, filename, 'start','end',0);

    %Plot_OCTCAP_local(path, filename, mouse, OCT, t_OCT, CAPs, t_CAP, Ampl, t_start, t_end, laser, genotype, saveyn)
end

%% - Make plot
clf
moving_avg_pts = 200;

subplot_tight(3,1,1)
plot(t_OCT,OCT{8,2}- movmean(OCT{8,2}, moving_avg_pts),'color',[0 0.5 0.9],'linewidth',2)
hold on
plot(t_OCT,OCT{6,2}- movmean(OCT{6,2}, moving_avg_pts),'k','linewidth',2)
text(19.8, 20, '50 dB')
text(19.8, 15, '70 dB','color',[0 0.5 0.9])
xlim([19.75 21.5])
ylim([-30 30])
title(filename,'interpreter','none')
ylabel('nm')

subplot_tight(3,1,2)
plot(t_OCT,OCT{1,2}- movmean(OCT{1,2}, moving_avg_pts),'k','linewidth',2)
text(29.8, 20, '34% laser')
xlim([29.75 31.5])
ylim([-30 30])

subplot_tight(3,1,3)
%plot(t_OCT,OCT{6,1} - movmean(OCT{6,1}, moving_avg_pts),'k','linewidth',2)
hold on
%plot(t_OCT,OCT{1,2}- movmean(OCT{1,2}, moving_avg_pts),'color',[0.8 0.8 0.8],'linewidth',3)
plot(t_OCT,OCT{6,2}- movmean(OCT{6,2}, moving_avg_pts),'color',[0.9 0.4 0.0],'linewidth',1.5)
plot(t_OCT,OCT{6,3}- movmean(OCT{6,3}, moving_avg_pts),'color',[0 0.5 0.9],'linewidth',1.5)
text(29.8, 20, '50 dB, 34% laser')
xlim([29.75 31.5])
ylim([-30 30])
xlabel('ms')

set(gcf,'position',[1000 0 450 1000])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_B.pdf', 'ContentType', 'vector');

%% - Fig Part C - subsites
% Set ylim for all plots here
ylimvec = [-40 40];
close all
% --- WT, 1/17, Mouse 1 --- %
% - Set the address of the mouse folder: 
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
exptdate = '011723';
mouse = 'ms1';

% - Enter the mouse folder
path = fullfile(folder,exptdate,mouse);
cd(path);

% - List all experiment files in the mouse folder
files = dir(path); files = {files.name}; files = files(4:end);

filesToPlot = 14; laser = 30;
filename = char(files(filesToPlot));
[t_OCT, OCT, ~, ~, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
offset_mean = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    subplot(5,1,1)
    plot(t_OCT,OCT{1,1}-offset_mean,'k','linewidth',2)
    line([laser laser], ylimvec,'color','r')
    xlim([laser-0.25 laser+1.5])
    ylim(ylimvec)
    set(gca, 'XTickLabel', [])
    ylabel('Bone')
    text(31,35, 'Bone')
    title('WT - 1/17/23 - Mouse1','interpreter','none')

filesToPlot = 15; laser = 30;
filename = char(files(filesToPlot));
[t_OCT, OCT, ~, ~, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
offset_mean = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    subplot(5,1,2)
    plot(t_OCT,OCT{1,1}-offset_mean,'k','linewidth',2)
    line([laser laser], ylimvec,'color','r')
    xlim([laser-0.25 laser+1.5])
    ylim(ylimvec)
    set(gca, 'XTickLabel', [])
    ylabel("Reissner's")
    text(31,35, "Reissner's")

filesToPlot = 13; laser = 30;
filename = char(files(filesToPlot));
[t_OCT, OCT, ~, ~, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
offset_mean = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    subplot(5,1,3)
    plot(t_OCT,OCT{1,1}-offset_mean,'k','linewidth',2)
    line([laser laser], ylimvec,'color','r')
    xlim([laser-0.25 laser+1.5])
    ylim(ylimvec)
    set(gca, 'XTickLabel', [])
    ylabel('BM')
    text(31,35, 'Basilar M.')

filesToPlot = 16; laser = 30;
filename = char(files(filesToPlot));
[t_OCT, OCT, ~, ~, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
offset_mean = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    subplot(5,1,4)
    plot(t_OCT,OCT{1,1}-offset_mean,'k','linewidth',2)
    line([laser laser], ylimvec,'color','r')
    xlim([laser-0.25 laser+1.5])
    ylim(ylimvec)
    set(gca, 'XTickLabel', [])
    ylabel('OHC')
    text(31,35, 'OHC')

filesToPlot = 58; laser = 40;
filename = char(files(filesToPlot));
[t_OCT, OCT, ~, ~, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
offset_mean = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    subplot(5,1,5)
    plot(t_OCT,OCT{1,1}-offset_mean,'k','linewidth',2)
    line([laser laser], ylimvec,'color','r')
    xlim([laser-0.25 laser+1.5])
    ylim(ylimvec)
    %set(gca, 'XTickLabel', [])
    ylabel('DEAD BM')
    text(41,35, 'Dead BM')
    xlabel('Time (ms)')

set(gcf,'Position',[1200 200 250 700])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_C_Subsites.pdf', 'ContentType', 'vector');

%% - Fig Part D - OCT & CAP overlaid
clear
%
%%%% --- ENTER REQUIRED INPUTS --- %%%%
% - Set the address of the mouse folder: 
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
exptdate = '011723';
mouse = 'ms2';

% - Enumerate the file #s that I wish to plot. Use the Excel file
% with experiment notes to decide which experiments to plot.
filesToPlot = 20;

% - Select plotting parameters
%t_start  = 29;      % Time in ms at which to start plotting the laser OCT and CAP
%t_end    = 34;      % Time in ms at which to stop plotting the laser OCT and CAP
%laser    = 30;      % Time at which laser was fired
laser    = 40;
t_start  = laser-1;
t_end    = laser + 4;
saveyn   = 0;       % 1 = save the plot in \Results\Oghalai\
genotype = 'WT';    % Enter mouse genotype (to properly label the plot and the file)

%%%% ------------------------------ %%%%

% - Enter the mouse folder
path = fullfile(folder,exptdate,mouse);
cd(path);

% - List all experiment files in the mouse folder
files = dir(path); files = {files.name}; files = files(4:end);

% - Loop over file to plot, importing their data and then plotting it
for f = 1 : length(filesToPlot)
    filename = char(files(filesToPlot(f)));

    [t_OCT, OCT, CAPs, t_CAP, Ampl] = SoundLaser_OCTCAP(path, filename, 'start','end',0);

    %Plot_OCTCAP(path, filename, mouse, OCT, t_OCT, CAPs, t_CAP, Ampl, t_start, t_end, laser, genotype, saveyn)
end

%v Plot
run = 6;
xlimvec = [59.5 64];

clf

subplot(3,1,1)
yyaxis left
plot(t_OCT, OCT{run,1} - 146.5,'k','linewidth',2)
ylim([-137 -113])
ylabel('nm')
line([60.1 60.1],[-137 -113],'color','k','linewidth',2)
line([61 61],[-137 -113],'color','k','linewidth',2)
grid on
title('2023_01_17 12_15_03 - Ms2 OCT & CAP','interpreter','none')

yyaxis right
plot(t_CAP+0.1, CAPs(run,:) -1,'r','linewidth',2)
xlim(xlimvec)
ylabel('uV')
ylim([-35 45])
grid on


subplot(3,1,2:3)
xlimvec = [60.1 61];
yyaxis left
plot(t_OCT, OCT{run,1} - 146.5,'k','linewidth',2)
ylim([-137 -113])
ylabel('nm')
grid on
%title('2023_01_17 12_15_03 - Ms2 OCT & CAP','interpreter','none')

yyaxis right
plot(t_CAP+0.1, CAPs(run,:) -1,'r','linewidth',2)
ylabel('uV')
xlim(xlimvec)
ylim([-35 45])
grid on

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_D_Overlaid_OCT_CAP.pdf', 'ContentType', 'vector');

%% - Fig part E, option 1 (1/17 Ms2) - Plot power vs CAP and OCT divot 
% --- WT, 1/17, Mouse 2 --- %%
% - Set the address of the mouse folder: 
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
exptdate = '011723';
mouse = 'ms2';

% - Enter the mouse folder
path = fullfile(folder,exptdate,mouse);
cd(path);

% - List all experiment files in the mouse folder
files = dir(path); files = {files.name}; files = files(4:end);

filesToPlot     = [37 38 39 40 41 42 43 44 45 46 47]; 
laser_on        = 30.*ones(1,12); % Time when laser comes on
laser_percent   = [34 10 13 16 15 14 17 18 20 26 34];

% Laser parameters
Target_Distance = 300;
Pulse_Length    = 100;
Fiber_Diameter  = 400;

% Calculate radiant exposure
[~,~,~,EnergyD] = CapellaEnergy(Fiber_Diameter, laser_percent, Target_Distance, Pulse_Length);

yLimVec = [-50 680];
offset_step = 60;

% Times within which to find OCT peak
OCT_min_start   = 30;
OCT_min_end     = 30.5;
OCT_max_start   = 30.35;
OCT_max_end     = 30.6;

% Define the times within which to find oCAP peak 1 amplitude
oCAP1_min_start   = 30; 
oCAP1_min_end     = 30.3;
oCAP1_max_start   = 30.2;
oCAP1_max_end     = 30.4;

% Define the times within which to find oCAP peak 2 amplitude
oCAP2_min_start   = 31; 
oCAP2_min_end     = 31.6;
oCAP2_max_start   = 31.5; 
oCAP2_max_end     = 32.5;

clf
offset = 0;

for  j = 1:length(filesToPlot)
filename = char(files(filesToPlot(j)));
[t_OCT, OCT, CAPs, t_CAP, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
subplot(1,2,1)
    offset_mean = mean(OCT{1,1}((laser_on-1)/0.01:laser_on/0.01));    
    signal = OCT{1,1}-offset_mean + offset;    
    plot(t_OCT, signal,'k','linewidth',2)
    hold on
        % Convert the times in ms to indices of the vector t
        [~, OCT_min_start_idx]  = findNearest(t_OCT, OCT_min_start);
        [~, OCT_min_end_idx]    = findNearest(t_OCT, OCT_min_end);
        [~, OCT_max_start_idx]  = findNearest(t_OCT, OCT_max_start);
        [~, OCT_max_end_idx]    = findNearest(t_OCT, OCT_max_end);

        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, OCT_min_start_idx, OCT_min_end_idx, OCT_max_start_idx, OCT_max_end_idx);
        OCT_Ampl(j) = Ampl;
        plot(t_OCT(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t_OCT(idx_max),signal(idx_max),'.g','markersize',20)
    xlim([laser_on(j)-0.5 laser_on(j)+5])
    ylim(yLimVec)
    title({filename},{'OCT'}, 'interpreter','none')
    ylabel('nm')
subplot(1,2,2)
    signal_CAP = CAPs + offset;
    plot(t_CAP,signal_CAP,'k','linewidth',2)
    hold on
        % Convert the times in ms to indices of the vector t
        [~, oCAP1_min_start_idx]  = findNearest(t_CAP, oCAP1_min_start);
        [~, oCAP1_min_end_idx]    = findNearest(t_CAP, oCAP1_min_end);
        [~, oCAP1_max_start_idx]  = findNearest(t_CAP, oCAP1_max_start);
        [~, oCAP1_max_end_idx]    = findNearest(t_CAP, oCAP1_max_end);
        
        % Convert the times in ms to indices of the vector t
        [~, oCAP2_min_start_idx]  = findNearest(t_CAP, oCAP2_min_start);
        [~, oCAP2_min_end_idx]    = findNearest(t_CAP, oCAP2_min_end);
        [~, oCAP2_max_start_idx]  = findNearest(t_CAP, oCAP2_max_start);
        [~, oCAP2_max_end_idx]    = findNearest(t_CAP, oCAP2_max_end);

        % Find peak 1
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal_CAP, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
        oCAP1_Ampl(j) = Ampl;
        % Markers for peak 1
        plot(t_CAP(idx_min),signal_CAP(idx_min),'.r','markersize',20)
        plot(t_CAP(idx_max),signal_CAP(idx_max),'.g','markersize',20)
        % Find peak 2
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal_CAP, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
        oCAP2_Ampl(j) = Ampl;
        % Markers for peak 2
        plot(t_CAP(idx_min),signal_CAP(idx_min),'.r','markersize',20)
        plot(t_CAP(idx_max),signal_CAP(idx_max),'.g','markersize',20)
    xlim([laser_on(j)-0.5 laser_on(j)+5])
    ylim(yLimVec)
    title('CAP')
    ylabel('uV')
    xlabel('ms')
offset = offset + offset_step;
    %subplot_tight(5,5,1,[0.02 0.02])
    %line([laser_on laser_on], ylimvec,'color','r')
end
set(gcf,'Position',[800 0 710 1000])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option1_Peaks_1_17_Ms2.pdf', 'ContentType', 'vector');


% - Plot OCT divot and oCAP1 divot vs Radiant energy
figure
yyaxis left
    plot(EnergyD,OCT_Ampl,'.','markersize',40);
    ylabel('OCT Divot Amplitude (nm)')
    title('OCT Divot vs oCAP1')
yyaxis right
    s = scatter(EnergyD,oCAP1_Ampl,'filled','SizeData',150);
    hold off
    alpha(s,0.7)
    ylabel('Quick divot amplitude (uV)')
    set(gca,'fontsize',16)
    xlabel('Radiant Energy (mJ/cm^2)')
set(gcf,'Position',[50 500 750 395])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option1_OCT_vs_oCAP1_1_17_Ms2.pdf', 'ContentType', 'vector');


% - Plot OCT divot, Action potential vs Radiant energy
figure
yyaxis left
    plot(EnergyD,OCT_Ampl,'.','markersize',40)
    ylabel('OCT Divot Amplitude (nm)')
    title('OCT Divot vs oCAP2')
yyaxis right
    plot(EnergyD,oCAP2_Ampl,'.','markersize',40)
    ylabel('CAP action potential amplitude (uV)')
    set(gca,'fontsize',16)
    xlabel('Radiant Energy (mJ/cm^2)')
set(gcf,'Position',[50 0 750 395])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option1_OCT_vs_oCAP2_1_17_Ms2.pdf', 'ContentType', 'vector');

% - oCAP peak 1 vs Action Potential
% clf
% 
% plot(oCAP1_Ampl, oCAP2_Ampl,'.k','markersize',40)
% xlabel('oCAP1 (uV)')
% ylabel('oCAP2 (uV)')
% set(gca, 'fontsize',16)
% title('CAP action potential (oCAP2) vs. Quick electrical divot (oCAP1) ')


%% - Fig part E, option 2 (1/17 Ms1) - Plot power vs CAP and OCT divot 
% --- WT, 1/17, Mouse 1 --- %%
% - Set the address of the mouse folder: 
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
exptdate = '011723';
mouse = 'ms1';

% - Enter the mouse folder
path = fullfile(folder,exptdate,mouse);
cd(path);

% - List all experiment files in the mouse folder
files = dir(path); files = {files.name}; files = files(4:end);

filesToPlot     = [26 27 28 29 30 33 34 37 39 40 41 42]; 
laser_on        = 40.*ones(1,12); % Time when laser comes on
laser_percent   = [34 10 13 16 20 18 16 17 27 51 68 100];

% Laser parameters
Target_Distance = 300;
Pulse_Length    = 100;
Fiber_Diameter  = 400;

% Calculate radiant exposure
[~,~,~,EnergyD] = CapellaEnergy(Fiber_Diameter, laser_percent, Target_Distance, Pulse_Length);

yLimVec = [-50 750];
offset_step = 60;

% Times within which to find OCT peak
OCT_min_start   = 40;
OCT_min_end     = 40.5;
OCT_max_start   = 40.35;
OCT_max_end     = 40.6;

% Define the times within which to find oCAP peak 1 amplitude
oCAP1_min_start   = 40; 
oCAP1_min_end     = 40.3;
oCAP1_max_start   = 40.2;
oCAP1_max_end     = 40.4;

% Define the times within which to find oCAP peak 2 amplitude
oCAP2_min_start   = 41; 
oCAP2_min_end     = 41.6;
oCAP2_max_start   = 41.5; 
oCAP2_max_end     = 42.5;

close all
offset = 0;

for  j = 1:length(filesToPlot)
filename = char(files(filesToPlot(j)));
[t_OCT, OCT, CAPs, t_CAP, ~] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
subplot(1,2,1)
    offset_mean = mean(OCT{1,1}((laser_on-1)/0.01:laser_on/0.01));    
    signal = OCT{1,1}-offset_mean + offset;    
    plot(t_OCT, signal,'k','linewidth',2)
    hold on
        % Convert the times in ms to indices of the vector t
        [~, OCT_min_start_idx]  = findNearest(t_OCT, OCT_min_start);
        [~, OCT_min_end_idx]    = findNearest(t_OCT, OCT_min_end);
        [~, OCT_max_start_idx]  = findNearest(t_OCT, OCT_max_start);
        [~, OCT_max_end_idx]    = findNearest(t_OCT, OCT_max_end);

        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, OCT_min_start_idx, OCT_min_end_idx, OCT_max_start_idx, OCT_max_end_idx);
        OCT_Ampl(j) = Ampl;
        plot(t_OCT(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t_OCT(idx_max),signal(idx_max),'.g','markersize',20)
    xlim([laser_on(j)-0.5 laser_on(j)+5])
    ylim(yLimVec)
    title({filename},{'OCT'}, 'interpreter','none')
    ylabel('nm')
subplot(1,2,2)
    signal_CAP = CAPs + offset;
    plot(t_CAP,signal_CAP,'k','linewidth',2)
    hold on
        % Convert the times in ms to indices of the vector t
        [~, oCAP1_min_start_idx]  = findNearest(t_CAP, oCAP1_min_start);
        [~, oCAP1_min_end_idx]    = findNearest(t_CAP, oCAP1_min_end);
        [~, oCAP1_max_start_idx]  = findNearest(t_CAP, oCAP1_max_start);
        [~, oCAP1_max_end_idx]    = findNearest(t_CAP, oCAP1_max_end);
        
        % Convert the times in ms to indices of the vector t
        [~, oCAP2_min_start_idx]  = findNearest(t_CAP, oCAP2_min_start);
        [~, oCAP2_min_end_idx]    = findNearest(t_CAP, oCAP2_min_end);
        [~, oCAP2_max_start_idx]  = findNearest(t_CAP, oCAP2_max_start);
        [~, oCAP2_max_end_idx]    = findNearest(t_CAP, oCAP2_max_end);

        % Find peak 1
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal_CAP, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
        oCAP1_Ampl(j) = Ampl;
        % Markers for peak 1
        plot(t_CAP(idx_min),signal_CAP(idx_min),'.r','markersize',20)
        plot(t_CAP(idx_max),signal_CAP(idx_max),'.g','markersize',20)
        % Find peak 2
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal_CAP, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
        oCAP2_Ampl(j) = Ampl;
        % Markers for peak 2
        plot(t_CAP(idx_min),signal_CAP(idx_min),'.r','markersize',20)
        plot(t_CAP(idx_max),signal_CAP(idx_max),'.g','markersize',20)
    xlim([laser_on(j)-0.5 laser_on(j)+5])
    ylim(yLimVec)
    title('CAP')
    ylabel('uV')
    xlabel('ms')
offset = offset + offset_step;
    %subplot_tight(5,5,1,[0.02 0.02])
    %line([laser_on laser_on], ylimvec,'color','r')
end
set(gcf,'Position',[800 0 710 1000])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option2_Peaks_1_17_Ms1.pdf', 'ContentType', 'vector');


% - Plot OCT divot and oCAP1 divot vs Radiant energy
figure
% fit_OCT = polyfit(EnergyD, OCT_Ampl,1);
% fit_CAP = polyfit(EnergyD, oCAP1_Ampl,1);
% x = linspace(0, max(EnergyD));

yyaxis left
    plot(EnergyD,OCT_Ampl,'.k','markersize',20);
    hold on
    %plot(x, x.*fit_OCT(1) + fit_OCT(2),'-k')
    ylabel('OCT Divot Amplitude (nm)')
    title('OCT Divot vs oCAP1')
    ylim([-20 200])
    set(gca,'Ycolor','k')
yyaxis right    
    plot(EnergyD,oCAP1_Ampl,'square','color','r','markersize',8,'linewidth',1.5);
    hold on
    %plot(x, x.*fit_CAP(1) + fit_CAP(2),'-r')
    %text(150, 20, num2str(fit_CAP(1)))
    hold off
    
    ylabel('oCAP1 divot amplitude (uV)')
    set(gca,'fontsize',16,'Ycolor','r')
    xlabel('Radiant Energy (mJ/cm^2)')
    ylim([-10 100])
set(gcf,'Position',[500 500 400 280])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option2_OCT_vs_oCAP1_1_17_Ms1.pdf', 'ContentType', 'vector');


% - Make it as two subplots that I cna later overlay
figure
subplot(2,1,1)
plot(EnergyD,OCT_Ampl,'.k','markersize',20);
    hold on
    %plot(x, x.*fit_OCT(1) + fit_OCT(2),'-k')
    ylabel('OCT Divot Amplitude (nm)')    
    ylim([-20 200])
    set(gca,'Ycolor','k','fontsize',16)

subplot(2,1,2)
    plot(EnergyD,oCAP1_Ampl,'square','color','r','markersize',8,'linewidth',1.5);          
    ylabel('oCAP1 divot amplitude (uV)')
    set(gca,'fontsize',16,'Ycolor','r')
    xlabel('Radiant Energy (mJ/cm^2)')
    ylim([-10 100])
set(gcf,'Position',[500 500 380 560])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option2_OCT_vs_oCAP1_1_17_Ms1_as_Subplots.pdf', 'ContentType', 'vector');



% - Plot OCT divot, Action potential vs Radiant energy
figure
yyaxis left
    plot(EnergyD,OCT_Ampl,'.','markersize',40)
    ylabel('OCT Divot Amplitude (nm)')
    title('OCT Divot vs oCAP2')
yyaxis right
    plot(EnergyD,oCAP2_Ampl,'.','markersize',40)
    ylabel('CAP action potential amplitude (uV)')
    set(gca,'fontsize',16)
    xlabel('Radiant Energy (mJ/cm^2)')
set(gcf,'Position',[50 0 750 395])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig1_E_Option2_OCT_vs_oCAP2_1_17_Ms1.pdf', 'ContentType', 'vector');

% - oCAP peak 1 vs Action Potential
% clf
% 
% plot(oCAP1_Ampl, oCAP2_Ampl,'.k','markersize',40)
% xlabel('oCAP1 (uV)')
% ylabel('oCAP2 (uV)')
% set(gca, 'fontsize',16)
% title('CAP action potential (oCAP2) vs. Quick electrical divot (oCAP1) ')

