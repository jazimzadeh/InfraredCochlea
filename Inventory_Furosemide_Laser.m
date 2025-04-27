%% - Make ABR timeline plot
% Plot thresholds, click wave 1, and laser against time (Mouse 1 2022_10_11)
clf
clear
exptdate = '2022_10_11';
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
path = strcat(folder,exptdate,'/');

% - oABRs 
cd(path)
clf

RunTwo          = readmatrix(strcat('Run2.csv'));
RunThree        = readmatrix(strcat('Run3.csv'));
RunFour         = readmatrix(strcat('Run4.csv'));
RunSix          = readmatrix(strcat('Run6.csv'));
RunEight        = readmatrix(strcat('Run8.csv'));
RunEleven       = readmatrix(strcat('Run11.csv'));
RunTwelve       = readmatrix(strcat('Run12.csv'));
RunFourteen     = readmatrix(strcat('Run14.csv'));
RunSixteen      = readmatrix(strcat('Run16.csv'));
RunSeventeen    = readmatrix(strcat('Run17.csv'));
RunNineteen     = readmatrix(strcat('Run19.csv'));
RunTwentyOne    = readmatrix(strcat('Run21.csv'));


Fs          = 25000;                    % Sampling rate in Hz
tstep       = 1/Fs;
t           = [0:(length(RunTwo)-49)]; % # of time points
t           = t.*tstep;                 % Time vector
laser_dur   = 0.150;                    % Laser pulse duration in ms

plot(t.*(1e3),(RunTwo(2,49:end)+8000),'linewidth',2)
hold on
plot(t.*(1e3),(RunThree(2,49:end)+6000),'linewidth',2)
plot(t.*(1e3),(RunFour(2,49:end)+4000),'linewidth',2)
plot(t.*(1e3),(RunSix(2,49:end)+2000),'linewidth',2)
plot(t.*(1e3),(RunEight(2,49:end)+0),'linewidth',2)
plot(t.*(1e3),(RunEleven(2,49:end)-2000),'linewidth',2)
plot(t.*(1e3),(RunTwelve(2,49:end)-4000),'linewidth',2)
plot(t.*(1e3),(RunFourteen(2,49:end)-6000),'linewidth',2)
plot(t.*(1e3),(RunSixteen(2,49:end)-8000),'linewidth',2)
plot(t.*(1e3),(RunSeventeen(2,49:end)-10000),'linewidth',2)
plot(t.*(1e3),(RunNineteen(2,49:end)-12000),'linewidth',2)
plot(t.*(1e3),(RunTwentyOne(2,49:end)-16000),'linewidth',2)

yLimMin = -18000;
yLimMax = 10000;
ylim([yLimMin yLimMax])
patch([1 1 (1+laser_dur) (1+laser_dur)],[yLimMin yLimMax yLimMax yLimMin],'r',...
    'EdgeColor', 'none', 'FaceAlpha', 0.3)

title(strcat(exptdate, ' -', ' Runs 2,3,4,6,8,11,12,14,16,17,19,21'), 'Interpreter', 'none')
xlabel('Time (ms)')
ylabel('Voltage nV')

legend('Run 2 - 1W, 150 us, 64 avg, 4 Hz', ...
    'Run 3 - 1W, 150 us, 64 avg, 4 Hz, changed pos.', ...
    'Run 4 - 2W, 150 us, 64 avg, 4 Hz', ...
    'Run 6 - 2W, 150 us, 64 avg, 4 Hz', ...
    'Run 8 - 5 min', ...
    'Run 11 - 16 min', ...
    'Run 12 - 19 min', ...
    'Run 14 - 24 min', ...
    'Run 16 - 33 min', ...
    'Run 17 - 35 min', ...
    'Run 19 - 41 min', ...
    'Run 21 - 49 min', ...
    'location','southeast')

% cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
% exportgraphics(gcf, strcat(exptdate,'_oABRs.pdf'), 'ContentType', 'vector');


% % % % % % % % % % % % % % % 
% - Manually get click ABR thresholds and timing w/respect to furosemide injection
click_times     = [-20 -8 2 7 12 20 27 37 45];
click_threshs   = [30 30 100 100 70 55 50 40 40];

laser_times     = [-10 -2 5 16 19 24 33 35 41 49];


% % % % % % % % % % % % % % % 
% - Timeline Plot: Click thresh, click 90 dB ABR, and oABRs on timeline

cd(path)
% Import the click runs
RunOne          = readmatrix(strcat('Run1.csv'));
RunFive         = readmatrix(strcat('Run5.csv'));
RunSeven        = readmatrix(strcat('Run7.csv'));
RunNine         = readmatrix(strcat('Run9.csv'));
RunTen          = readmatrix(strcat('Run10.csv'));
RunThirteen     = readmatrix(strcat('Run13.csv'));
RunFifteen      = readmatrix(strcat('Run15.csv'));
RunEighteen     = readmatrix(strcat('Run18.csv'));
RunTwenty       = readmatrix(strcat('Run20.csv'));

clf
subplot(1,6,1:2)
plot(click_threshs, -click_times,'.k','markersize',20)
ylim([-51 25])
xlim([0 100])
xlabel('Threshold (dB SPL)')
ylabel('Time (min)')
grid on
title('10/11/2022 Mouse 1')

subplot(1,6,3:4)
plot(t.*(1e3),(RunOne(2,49:end)+20000),'linewidth',2,'color',[0 0 0])
hold on
plot(t.*(1e3),(RunFive(2,49:end)+8000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunSeven(2,49:end)-2000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunNine(2,49:end) -7000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunTen(2,49:end)  -12000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunThirteen(2,49:end)-20000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunFifteen(2,49:end) -27000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunEighteen(2,49:end)-37000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunTwenty(2,49:end)  -45000),'linewidth',2,'color',[0 0 0])
yLimMin = -51000;
yLimMax = 25000;
ylim([yLimMin yLimMax])
title('90 dB Click ABR')

% The y-offset is the time (min) x 1000. Y scale should be like subplot 1
% but x 1000.
subplot(1,6,5:6)
plot(t.*(1e3),(RunFour(2,49:end)+10000),'linewidth',2,'color',[0 0 0])
hold on
plot(t.*(1e3),(RunSix(2,49:end)+2000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunEight(2,49:end)-5000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunEleven(2,49:end)-16000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunTwelve(2,49:end)-19000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunFourteen(2,49:end)-24000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunSeventeen(2,49:end)-35000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunNineteen(2,49:end)-41000),'linewidth',2,'color',[0 0 0])
plot(t.*(1e3),(RunTwentyOne(2,49:end)-49000),'linewidth',2,'color',[0 0 0])
title('Laser ABR')

ylim([yLimMin yLimMax])
set(gcf,'Position',[1000 800 750 300])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Compiled_Furosemide_ABR_Timeline.pdf', 'ContentType', 'vector');

%% - Make structure w/data from all mice
clear
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets")
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";

Table       = readtable('Inventory_Furosemide_Laser.csv');
TableHeader = Table.Properties.VariableNames;

for j = 1:size(Table,1)   % Iterate through the table rows (mice); create Data structure
    Data(j).Date            = Table.Date(j);
    Data(j).Mouse           = Table.Mouse(j);
    Data(j).Strain          = Table.Strain(j);

    DataPath                = cd(strcat(folder,Table.Date(j)));    
    Data(j).Pre_click       = readmatrix(strcat('Run',num2str(Table.Pre_click(j)),'.csv'));
    Data(j).Pre_oCAP        = readmatrix(strcat('Run',num2str(Table.Pre_oCAP(j)),'.csv'));
    Data(j).During_click    = readmatrix(strcat('Run',num2str(Table.During_click(j)),'.csv'));
    Data(j).During_oCAP     = readmatrix(strcat('Run',num2str(Table.During_oCAP(j)),'.csv'));
    Data(j).Post_click      = readmatrix(strcat('Run',num2str(Table.Post_click(j)),'.csv'));
    Data(j).Post_oCAP       = readmatrix(strcat('Run',num2str(Table.Post_oCAP(j)),'.csv'));
end

% Save the Data structure to a mat file
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
save("Furosemide_Laser_Inventory.mat","Data")

%% - Load the Data structure 
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
load ("Furosemide_Laser_Inventory.mat")

t  = linspace(0,0.01098,268).*1000;

%% - Plot click and laser response: Pre, During, and Post Furosemide
clf

offset_step = 10;
yLimVec = [-5 25];
xLimVec = [0 10];

% Define the times within which to find wave 1 amplitude
min_start   = 2; % Input time (in ms) at which to start looking for minimum
min_end     = 3;
max_start   = 1.8;
max_end     = 2.5;

% Convert the times in ms to indices of the vector t
[~, min_start_idx]  = findNearest(t, min_start);
[~, min_end_idx]    = findNearest(t, min_end);
[~, max_start_idx]  = findNearest(t, max_start);
[~, max_end_idx]    = findNearest(t, max_end);

subplot_tight(1,6,1)
offset = 0;
for p = 1:size(Data,2)
    signal = 0.001.*Data(p).Pre_click(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on    
    text(4.5,offset+3, strcat(Data(p).Date, ' - Ms ', num2str(Data(p).Mouse)),'interpreter','none')
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    Pre_Click_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_step;
    grid on
    ylabel('uV')
    xlabel('ms')
    ylim(yLimVec)
    xlim(xLimVec)
end
title('Click, Pre, 90 dB')

subplot_tight(1,6,2)
offset = 0;
for p = 1:size(Data,2)
    signal = 0.001.*Data(p).Pre_oCAP(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on        
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    Pre_Laser_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_step;
    grid on
    xlabel('ms')
    ylim(yLimVec)
    xlim(xLimVec)
end
title('Laser, Pre, 34%')

subplot_tight(1,6,3)
offset = 0;
for p = 1:size(Data,2)
    signal = 0.001.*Data(p).During_click(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on        
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    During_Click_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_step;
    grid on
    xlabel('ms')
    ylim(yLimVec)
    xlim(xLimVec)
end
title('Click, During, 90 dB')

subplot_tight(1,6,4)
offset = 0;
for p = 1:size(Data,2)
    signal = 0.001.*Data(p).During_oCAP(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on        
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    During_Laser_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_step;
    grid on
    xlabel('ms')
    ylim(yLimVec)
    xlim(xLimVec)
end
title('Laser, During, 34%')

subplot_tight(1,6,5)
offset = 0;
for p = 1:size(Data,2)
    signal = 0.001.*Data(p).Post_click(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on        
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    Post_Click_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_step;
    grid on
    xlabel('ms')
    ylim(yLimVec)
    xlim(xLimVec)
end
title('Click, Post, 90 dB')

subplot_tight(1,6,6)
offset = 0;
for p = 1:size(Data,2)
    signal = 0.001.*Data(p).Post_oCAP(2,49:316) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on        
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, min_start_idx, min_end_idx, max_start_idx, max_end_idx);
    Post_Laser_Ampl(p) = Ampl;
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    offset = offset + offset_step;
    grid on
    xlabel('ms')
    ylim(yLimVec)
    xlim(xLimVec)
end
title('Laser, Post, 34%')

set(gcf,'Position',[0 500 1500 400])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_Furosemide_Laser_PreDuringPost.pdf'), 'ContentType', 'vector');

%% - Histograms of Wave 1 Amplitudes (Pre, During, Post Furosemide) for click and laser stim
clf
subplot(1,2,1)
plot(0.7+0.4.*(rand(1,length(Pre_Click_Ampl))),Pre_Click_Ampl,'.k','markersize',40)
hold on
plot(2.7+0.4.*(rand(1,length(During_Click_Ampl))),During_Click_Ampl,'.k','markersize',40)
plot(4.7+0.4.*(rand(1,length(Post_Click_Ampl))),Post_Click_Ampl,'.k','markersize',40)
xlim([0 6])
ylim([0 10])
xticks([1 3 5])
set(gca,'fontsize',16,'xticklabel',{'Pre', 'During furosemide', 'Post'})
title('Click 90 dB - Wave 1 Ampl')
ylabel('uV')
%[H, p] = paired_ttest(Pre_Click_Ampl, During_Click_Ampl, 1, 9, [1 3]);
%[H, p] = paired_ttest(Pre_Click_Ampl, Post_Click_Ampl, 1, 9.7, [1 5]);
%[H, p] = paired_ttest(During_Click_Ampl, Post_Click_Ampl, 1, 10.5, [3 5]);

[p] = KruskalWallisPlot(Pre_Click_Ampl', During_Click_Ampl', 1, 9, [1 3]);
[p] = KruskalWallisPlot(Pre_Click_Ampl', Post_Click_Ampl', 1, 9.7, [1 5]);
[p] = KruskalWallisPlot(During_Click_Ampl', Post_Click_Ampl', 1, 3, [3 5]);
% [H, p] = ttest(Pre_Click_Ampl, During_Click_Ampl);
% line([1 3],[9 9],'linewidth',2,'color','k')
% if H == 1
%     if p < 0.001
%         text(1.9,9.1,'***','fontsize',40)
%     elseif p < 0.01
%         text(1.9,9.1,'**','fontsize',40)
%     elseif p < 0.05
%         text(1.9,9.1,'*','fontsize',40)
%     end
% else
%     text(2,9.1,'NS')
% end

subplot(1,2,2)
plot(0.7+0.4.*(rand(1,length(Pre_Laser_Ampl))),Pre_Laser_Ampl,'.k','markersize',40)
hold on
plot(2.7+0.4.*(rand(1,length(During_Laser_Ampl))),During_Laser_Ampl,'.k','markersize',40)
plot(4.7+0.4.*(rand(1,length(Post_Laser_Ampl))),Post_Laser_Ampl,'.k','markersize',40)
xlim([0 6])
ylim([0 10])
xticks([1 3 5])
set(gca,'fontsize',16,'xticklabel',{'Pre', 'During furosemide', 'Post'})
title('Laser 34% - Wave 1 Ampl')
ylabel('uV')
%[H, p] = paired_ttest(Pre_Laser_Ampl, Post_Laser_Ampl, 1, 8, [1 5]);
%[H, p] = paired_ttest(Pre_Laser_Ampl, During_Laser_Ampl, 1, 9, [1 3]);
%[H, p] = paired_ttest(Post_Laser_Ampl, During_Laser_Ampl, 1, 7, [3 5]);

[p] = KruskalWallisPlot(Pre_Laser_Ampl', Post_Laser_Ampl', 1, 8, [1 5]);
[p] = KruskalWallisPlot(Pre_Laser_Ampl', During_Laser_Ampl', 1, 9, [1 3]);
[p] = KruskalWallisPlot(Post_Laser_Ampl', During_Laser_Ampl', 1, 7, [3 5]);

set(gcf,'Position', [1150 800 550 300])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Compiled_Furosemide_Laser_Histogram.pdf'), 'ContentType', 'vector');

%% - Save the amplitudes into a CSV table 
Pre_Click_Ampl      = Pre_Click_Ampl';
During_Click_Ampl   = During_Click_Ampl';
Post_Click_Ampl     = Post_Click_Ampl';
Pre_Laser_Ampl      = Pre_Laser_Ampl';
During_Laser_Ampl   = During_Laser_Ampl';
Post_Laser_Ampl     = Post_Laser_Ampl';

Summary = table(Pre_Click_Ampl, During_Click_Ampl, Post_Click_Ampl, Pre_Laser_Ampl, During_Laser_Ampl, Post_Laser_Ampl);

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
writetable(Summary,'Furosemide_Laser_Wave1.csv')

%% Two-sample paired t-test
% I confirmed this with Rstudio, it gives me the same p-values
clc
[H, p] = ttest(Pre_Click_Ampl, During_Click_Ampl)
[H, p] = ttest(Pre_Laser_Ampl, During_Laser_Ampl)

