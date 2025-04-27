%% - Part A - Sound and Laser, OCT & CAP: WT, Tecta, VGLUT3 time traces

close all
clear

offset_CAP = 60;
yLimOCT         = [-70 70];
yLimCAP         = [-55 55];

subplot(9,2,1)
    %%%% --- WT, 1/15, Mouse 1 --- %%%%
    % - Set the address of the mouse folder: 
    folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
    exptdate = '011523';
    mouse = 'ms1';
    
    % - Enter the mouse folder
    path = fullfile(folder,exptdate,mouse);
    cd(path);
    
    % - List all experiment files in the mouse folder
    files = dir(path); files = {files.name}; files = files(4:end);
    
    filesToPlot = 11; 
    laser = 30; % Laser onset
    sound = 20; % Sound onset
    filename = char(files(filesToPlot));
    [t_OCT, OCT, CAPs, t_CAP, Amplitudes] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
    idx_70dB = find(Amplitudes == 70);
    idx_40dB = find(Amplitudes == 40);
    offset_sound70  = mean(OCT{idx_70dB,2}((sound-1)/0.01:sound/0.01));
    offset_sound40  = mean(OCT{idx_40dB,2}((sound-1)/0.01:sound/0.01));
    offset_laser    = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    
    plot(t_OCT,OCT{idx_70dB,2}-offset_sound70,'k','linewidth',1.5)
    text(sound-0.4,25,'70dB')
    
    xlim([sound-0.5 sound+2])
    ylim(yLimOCT) 
    ylabel('nm')
    title({'OCT',filename},'interpreter','none')

subplot(9,2,3)
    plot(t_OCT,OCT{idx_40dB,2}-offset_sound40,'k','linewidth',1.5)
        text(sound-0.4,25,'40dB')
        xlim([sound-0.5 sound+2])
        ylim(yLimOCT) 
        ylabel('WT')

subplot(9,2,5)
    plot(t_OCT,OCT{1,1}-offset_laser,'k','linewidth',1.5)
    text(laser-0.4,25,'Laser')
    xlim([laser-0.5 laser+2])
    ylim(yLimOCT)

subplot(9,2,2) % - CAPs
    plot(t_CAP+0.1, CAPs(idx_70dB,:),'k','linewidth',1.5)       
    xlim([sound-1 sound+9])
    ylim(yLimCAP)
    ylabel('uV')
    title('CAPs')

subplot(9,2,4)
    plot(t_CAP+0.1, CAPs(idx_40dB,:),'k','linewidth',1.5)
        xlim([sound-1 sound+9])
        ylim(yLimCAP)

subplot(9,2,6)
    plot(t_CAP+0.1, CAPs(1,:),'k','linewidth',1.5)
    xlim([laser-1 laser+9])
    ylim(yLimCAP)



%%%% --- Tecta, 1/18, Mouse 1 --- %%%%
subplot(9,2,7)    
    % - Set the address of the mouse folder: 
    folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
    exptdate = '011823';
    mouse = 'ms2';
    
    % - Enter the mouse folder
    path = fullfile(folder,exptdate,mouse);
    cd(path);
    
    % - List all experiment files in the mouse folder
    files = dir(path); files = {files.name}; files = files(4:end);
    
    filesToPlot = 10; 
    laser = 30; % Laser onset
    sound = 20; % Sound onset
    filename = char(files(filesToPlot));
    [t_OCT, OCT, CAPs, t_CAP, Amplitudes] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
    idx_70dB = find(Amplitudes == 70);
    idx_40dB = find(Amplitudes == 40);
    offset_sound70  = mean(OCT{idx_70dB,2}((sound-1)/0.01:sound/0.01));
    offset_sound40  = mean(OCT{idx_40dB,2}((sound-1)/0.01:sound/0.01));
    offset_laser    = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    
    plot(t_OCT,OCT{idx_70dB,2}-offset_sound70,'k','linewidth',1.5)
    text(sound-0.4,25,'70dB')
    
    xlim([sound-0.5 sound+2])
    ylim(yLimOCT) 
    ylabel('nm')
    title(filename,'interpreter','none')

subplot(9,2,9)
    plot(t_OCT,OCT{idx_40dB,2}-offset_sound40,'k','linewidth',1.5)
        text(sound-0.4,25,'40dB')
        xlim([sound-0.5 sound+2])
        ylim(yLimOCT) 
        ylabel('Tecta')

subplot(9,2,11)
    plot(t_OCT,OCT{1,1}-offset_laser,'k','linewidth',1.5)
    text(laser-0.4,25,'Laser')
    xlim([laser-0.5 laser+2])
    ylim(yLimOCT)

subplot(9,2,8) % - CAPs
    plot(t_CAP+0.1, CAPs(idx_70dB,:),'k','linewidth',1.5)       
    xlim([sound-1 sound+9])
    ylim(yLimCAP)
    ylabel('uV')

subplot(9,2,10)
    plot(t_CAP+0.1, CAPs(idx_40dB,:),'k','linewidth',1.5)
        xlim([sound-1 sound+9])
        ylim(yLimCAP)

subplot(9,2,12)
    plot(t_CAP+0.1, CAPs(1,:),'k','linewidth',1.5)
    xlim([laser-1 laser+9])
    ylim(yLimCAP)


%%%% --- VGLUT3, 1/18, Mouse 4 --- %%%%
subplot(9,2,13)    
    % - Set the address of the mouse folder: 
    folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";
    exptdate = '011823';
    mouse = 'ms4';
    
    % - Enter the mouse folder
    path = fullfile(folder,exptdate,mouse);
    cd(path);
    
    % - List all experiment files in the mouse folder
    files = dir(path); files = {files.name}; files = files(4:end);
    
    filesToPlot = 8; 
    laser = 30; % Laser onset
    sound = 20; % Sound onset
    filename = char(files(filesToPlot));
    [t_OCT, OCT, CAPs, t_CAP, Amplitudes] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
    idx_70dB = find(Amplitudes == 70);
    idx_40dB = find(Amplitudes == 40);
    offset_sound70  = mean(OCT{idx_70dB,2}((sound-1)/0.01:sound/0.01));
    offset_sound40  = mean(OCT{idx_40dB,2}((sound-1)/0.01:sound/0.01));
    offset_laser    = mean(OCT{1,1}((laser-1)/0.01:laser/0.01));
    
    plot(t_OCT,OCT{idx_70dB,2}-offset_sound70,'k','linewidth',1.5)
    text(sound-0.4,25,'70dB')
    
    xlim([sound-0.5 sound+2])
    ylim(yLimOCT) 
    ylabel('nm')
    title(filename,'interpreter','none')

subplot(9,2,15)
    plot(t_OCT,OCT{idx_40dB,2}-offset_sound40,'k','linewidth',1.5)
        text(sound-0.4,25,'40dB')
        xlim([sound-0.5 sound+2])
        ylim(yLimOCT) 
        ylabel('VGLUT3')

subplot(9,2,17)
    plot(t_OCT,OCT{1,1}-offset_laser,'k','linewidth',1.5)
    text(laser-0.4,25,'Laser')
    xlim([laser-0.5 laser+2])
    ylim(yLimOCT)
    xlabel('ms')

subplot(9,2,14) % - CAPs
    % Filter at 10 kHz (we had changed the filter to 50 kHz during expt),
    % others were 10 kHz
    Hd = designfilt('lowpassfir','FilterOrder',20,'CutoffFrequency',10000, ...
       'DesignMethod','window','Window',{@kaiser,3},'SampleRate',200000);
    y1 = filter(Hd,CAPs(idx_70dB,:));
    %plot(t_CAP+0.1, CAPs(idx_70dB,:),'k','linewidth',1.5)  
    plot(t_CAP+0.1, y1,'k','linewidth',1.5)  
    xlim([sound-1 sound+9])
    ylim(yLimCAP)
    ylabel('uV')

subplot(9,2,16)
    y1 = filter(Hd,CAPs(idx_40dB,:));
    plot(t_CAP+0.1, y1,'k','linewidth',1.5)
    %plot(t_CAP+0.1, CAPs(idx_40dB,:),'k','linewidth',1.5)
        xlim([sound-1 sound+9])
        ylim(yLimCAP)

subplot(9,2,18)
    Hd = designfilt('lowpassfir','FilterOrder',20,'CutoffFrequency',10000, ...
       'DesignMethod','window','Window',{@kaiser,3},'SampleRate',200000);
    y1 = filter(Hd,CAPs(1,:));
    %plot(t_CAP+0.1, CAPs(1,:),'k','linewidth',1.5)
    plot(t_CAP+0.1, y1,'k','linewidth',1.5)
    xlim([laser-1 laser+9])   
    ylim(yLimCAP)
    xlabel('ms')

set(gcf,'Position',[1000 0 500 1200])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figures_OCT_Fig2_A_Mutant_OCT_CAPs.pdf', 'ContentType', 'vector');

%% - Part B - Make structure w/data from all mice (all 12 mices from all 3 genotypes)
clear
close all
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets")
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";

Table       = readtable('Inventory_OCT_Mutants_Laser_Motion.csv');
TableHeader = Table.Properties.VariableNames;
folder      = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";

for j = 1:size(Table,1)   % Iterate through the table rows (mice); create Data structure
    Data(j).Date            = Table.exptdate(j);
    Data(j).Mouse           = Table.mouse(j);
    Data(j).Genotype_WT0_Tecta1_VGLUT2= Table.genotype_WT0_Tecta1_VGLUT2(j);
    path                    = fullfile(folder, strcat('0',num2str(Table.exptdate(j))), Table.mouse(j));
    files = dir(path); files = {files.name}; files = files(4:end);
    filename    = char(files(Table.filenumber(j)));
    file_tuning = char(files(Table.file_tuning(j)));
    [Data(j).t_OCT, Data(j).OCT, Data(j).CAPs, Data(j).t_CAP, Data(j).Amplitudes] = SoundLaser_OCTCAP(path, filename, 'start','end',0);
    Data(j).laser_onset     = Table.laser_onset(j);
    Data(j).laser_power     = Table.laser_power(j);
    Data(j).filename        = filename;
    Data(j).file_tuning     = file_tuning;
end

% Save the Data structure to a mat file
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
save("OCT_Mutants_Inventory.mat","Data")


%% - Load the Data structure with OCT traces for each mutant type (WT, Tecta, VGLUT3)
clear
close all
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Inventory_Workspace_Data")
load ("OCT_Mutants_Inventory.mat")

Data_WT     = Data([Data.Genotype_WT0_Tecta1_VGLUT2]==0);
Data_Tecta  = Data([Data.Genotype_WT0_Tecta1_VGLUT2]==1);
Data_VGLUT  = Data([Data.Genotype_WT0_Tecta1_VGLUT2]==2);

%% - Show effect of baseline subtraction (red) and moving avg subtraction (black)
clf

for p = 1:size(Data,2)
    subplot_tight(6,2,p)
        signal      = Data(p).OCT{1,1};
        baseline    = mean(signal((Data(p).laser_onset-0.5)/.01:(Data(p).laser_onset/0.01)));
        moving_avg  = movmean(signal,200);
        % Original:
        %plot(Data(p).t_OCT, signal,'color','k','linewidth',2) 
        % Plot after baseline subtraction:
        plot(Data(p).t_OCT, signal-baseline,'color','r','linewidth',2)
        hold on
        plot(Data(p).t_OCT, signal-moving_avg,'color','k','linewidth',2)
        text((Data(p).laser_onset - 1.85), 55, genotype_converter(Data(p).Genotype_WT0_Tecta1_VGLUT2), 'fontsize', 18)
        xlim([(Data(p).laser_onset - 2) (Data(p).laser_onset +5)])
        ylim([-70 75])

        if p == 1
            title('Red = Baseline-subtracted, Black = Moving avg subtracted')
        end

        if p == 2
            title('Y = Displacement (nm), X = Time (ms)')
        end        
end

set(gcf,'Position',[700 0 800 1000])


%% - Calculate V_RMS and overlay it on the moving avg corrected plots
clf
clear RMS_Tecta RMS_VGLUT3 RMS_WT

for p = 1:size(Data,2)
    subplot_tight(6,2,p)
        signal      = Data(p).OCT{1,1};       
        moving_avg  = movmean(signal,200);
        signal_z    = signal - moving_avg;
        plot(Data(p).t_OCT, signal_z,'color','k','linewidth',2)
        hold on
        
        start_RMS   = Data(p).laser_onset;
        stop_RMS    = Data(p).laser_onset + 2;
        signal_toRMS= signal_z((start_RMS/0.01) : (stop_RMS/0.01));
        signal_RMS(p) = sqrt(sum(signal_toRMS.^2));  

        text((Data(p).laser_onset - 0.85), 50, genotype_converter(Data(p).Genotype_WT0_Tecta1_VGLUT2), 'fontsize', 16)
        text((Data(p).laser_onset + 2), 50, strcat('V_{RMS} = ',num2str(round(signal_RMS(p))),' uV'), 'fontsize', 16)
        text((Data(p).laser_onset - 0.85), -55, Data(p).filename,'interpreter','none')
        xlim([(Data(p).laser_onset - 1) (Data(p).laser_onset +4)])
        ylim([-70 70])

        if Data(p).Genotype_WT0_Tecta1_VGLUT2 == 0
            RMS_WT(p) = signal_RMS(p);
        elseif Data(p).Genotype_WT0_Tecta1_VGLUT2 == 1
            RMS_Tecta(p) = signal_RMS(p);
        elseif Data(p).Genotype_WT0_Tecta1_VGLUT2 == 2
            RMS_VGLUT3(p) = signal_RMS(p);
        end

        if p == 1
            title('Red = Baseline-subtracted, Black = Moving avg subtracted')
        end

        if p == 2
            title('Y = Displacement (nm), X = Time (ms)')
        end        
end

RMS_WT      = nonzeros(RMS_WT);
RMS_Tecta   = nonzeros(RMS_Tecta);
RMS_VGLUT3  = nonzeros(RMS_VGLUT3);

set(gcf,'Position',[700 0 800 1000])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_B_Mutant_OCT_RMS_TimeSeries.pdf'), 'ContentType', 'vector');

%% - Histogram of RMS OCT values of WT, Tecta, VGLUT3
clf
plot(0.7+0.4.*(rand(1,length(RMS_WT))), RMS_WT,'.k','markersize',40)
hold on
plot(2.7+0.4.*(rand(1,length(RMS_Tecta))), RMS_Tecta,'.k','markersize',40)
plot(4.7+0.4.*(rand(1,length(RMS_VGLUT3))), RMS_VGLUT3,'.k','markersize',40)

%[H, p] = unpaired_ttest(RMS_WT, RMS_Tecta, 1, 260, [1 3]);
%[H2, p2] = unpaired_ttest(RMS_WT, RMS_VGLUT3, 1, 285, [1 5]);
%[H3, p3] = unpaired_ttest(RMS_Tecta, RMS_VGLUT3, 1, 273, [3 5]);

[p] = KruskalWallisPlot(RMS_WT, RMS_Tecta, 1, 260, [1 3]);
[p] = KruskalWallisPlot(RMS_WT, RMS_VGLUT3, 1, 285, [1 5]);
[p] = KruskalWallisPlot(RMS_Tecta, RMS_VGLUT3, 1, 273, [3 5]);

xticks([1 3 5])
xlim([0 6])
set(gca,'fontsize',16,'xticklabel',{'WT', 'Tecta', 'VGLUT3'})
ylabel('RMS (uV)')

set(gcf,'Position',[700 300 600 500])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_B_Mutant_OCT_RMS_Histogram.pdf'), 'ContentType', 'vector');

%% - Plot trace, spectrogram, power spectrum -> Wild Type
clf
for p = 1:size(Data_WT,2)
    subplot_tight(7,size(Data_WT,2), p)
        signal      = Data_WT(p).OCT{1,1};    
        idx_baseline_start  = (Data_WT(p).laser_onset - 0.5)/0.01;
        idx_baseline_end    = (Data_WT(p).laser_onset)/0.01;
        baseline    = mean(signal(idx_baseline_start : idx_baseline_end));
        moving_avg  = movmean(signal,200);        
        signal_z    = signal - moving_avg;
        signal_b    = signal - baseline;

        plot(Data_WT(p).t_OCT, signal_z,'color','k','linewidth',2)
        text((Data_WT(p).laser_onset - 0.85), 55, Data_WT(p).filename,'interpreter','none','fontsize',8)
        title('WT')
        xlim([(Data_WT(p).laser_onset - 1) (Data_WT(p).laser_onset +3)])
        ylim([-70 70])

    subplot_tight(7, size(Data_WT,2), [(p+5) (p+15)])
        window   = 70;      % 50
        noverlap = window-1;      % 75% of window length; 30
        nfft     = 10000;   % # of frequencies
        Fs       = 100000;  % Sampling rate

        spectrogram(signal_b ,window ,noverlap ,nfft ,Fs ,'yaxis'); % Short-time fourier transform
        colormap(flipud(turbo))
        clim([-45 0])
        ylim([0 15])
        xlim([(Data_WT(p).laser_onset - 1) (Data_WT(p).laser_onset +2)])

    subplot_tight(7, size(Data_WT,2), [(p+20) (p+30)])
        t_start     = Data_WT(p).laser_onset; % in ms
        t_end       = (Data_WT(p).laser_onset) + 2;
        Fs          = 100000;
        [power_spec, freqs] = pspectrum(Data_WT(p).OCT{1,1}(t_start/0.01 : t_end/0.01), Fs, 'power');
        plot(freqs, power_spec, 'k','linewidth',2)
        ylim([0 1])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        hold on

        [troughs, troughlocs]   = findpeaks(-power_spec); 
        idx_first_trough        = troughlocs(1);
        [pk_max, idx_pk_max]    = findpeaks(power_spec, 'SortStr', 'descend','NPeaks',1);
        [~ , idx_high_troughs]  = findpeaks(-power_spec(idx_pk_max:end));
        idx_high_trough         = idx_pk_max + idx_high_troughs(1);
        [~ , idx_low_troughs]   = findpeaks(-power_spec(1:idx_pk_max));
        idx_low_trough          = idx_low_troughs(end);
        pow = sum(power_spec(idx_low_trough:idx_high_trough)) * freqs(2);

        plot(freqs(1:idx_first_trough), power_spec(1:idx_first_trough), 'r','linewidth',2)    
        plot(freqs(idx_pk_max), power_spec(idx_pk_max),'.','color',[0 0.8 0.1],'markersize',30)           
        plot(freqs(idx_high_trough), power_spec(idx_high_trough), '.','color',[0 0.4 1],'markersize',30)       
        plot(freqs(idx_low_trough), power_spec(idx_low_trough), '.','color',[0 0.4 1],'markersize',30)       
        text(10000, 0.8, strcat(num2str(round(pow)), ' uV^2'),'fontsize',14)
        
        Powers_WT(p)    = pow;
        Peak_freq_WT(p) = freqs(idx_pk_max); 

        ylim([0 1])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        if p == 1
            ylabel('Power Density $(\mu V^2 / Hz)$','Interpreter','latex','fontsize',14)
        end    
        
end

set(gcf,'Position',[0 100 1500 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_WT_Spectrograms_PowerSpecs.pdf'), 'ContentType', 'vector');

%% - Plot trace, spectrogram, power spectrum -> Tecta
clf
for p = 1:size(Data_Tecta,2)
    subplot_tight(7,size(Data_Tecta,2), p)
        signal      = Data_Tecta(p).OCT{1,1};    
        idx_baseline_start  = (Data_Tecta(p).laser_onset - 0.5)/0.01;
        idx_baseline_end    = (Data_Tecta(p).laser_onset)/0.01;
        baseline    = mean(signal(idx_baseline_start : idx_baseline_end));
        moving_avg  = movmean(signal,200);        
        signal_z    = signal - moving_avg;
        signal_b    = signal - baseline;

        plot(Data_Tecta(p).t_OCT, signal_z,'color','k','linewidth',2)
        text((Data_Tecta(p).laser_onset - 0.85), 55, Data_Tecta(p).filename,'interpreter','none','fontsize',8)
        title('Tecta')
        xlim([(Data_Tecta(p).laser_onset - 1) (Data_Tecta(p).laser_onset +3)])
        ylim([-70 70])

    subplot_tight(7, size(Data_Tecta,2), [(p+4) (p+12)])
        window   = 70;      % 50
        noverlap = window-1;      % 75% of window length; 30
        nfft     = 10000;   % # of frequencies
        Fs       = 100000;  % Sampling rate

        spectrogram(signal_b ,window ,noverlap ,nfft ,Fs ,'yaxis'); % Short-time fourier transform
        colormap(flipud(turbo))
        clim([-45 0])
        ylim([0 15])
        xlim([(Data_Tecta(p).laser_onset - 1) (Data_Tecta(p).laser_onset +2)])

    subplot_tight(7, size(Data_Tecta,2), [(p+16) (p+24)])
        t_start     = Data_Tecta(p).laser_onset; % in ms
        t_end       = (Data_Tecta(p).laser_onset) + 2;
        Fs          = 100000;
        [power_spec, freqs] = pspectrum(Data_Tecta(p).OCT{1,1}(t_start/0.01 : t_end/0.01), Fs);
        plot(freqs, power_spec, 'k','linewidth',2)
        hold on

        [troughs, troughlocs]   = findpeaks(-power_spec); 
        idx_first_trough        = troughlocs(1);
        [pk_max, idx_pk_max]    = findpeaks(power_spec, 'SortStr', 'descend','NPeaks',1);
        [~ , idx_high_troughs]  = findpeaks(-power_spec(idx_pk_max:end));
        idx_high_trough         = idx_pk_max + idx_high_troughs(1);
        [~ , idx_low_troughs]   = findpeaks(-power_spec(1:idx_pk_max));
        idx_low_trough          = idx_low_troughs(end);
        pow = sum(power_spec(idx_low_trough:idx_high_trough)) * freqs(2);

        plot(freqs(1:idx_first_trough), power_spec(1:idx_first_trough), 'r','linewidth',2)    
        plot(freqs(idx_pk_max), power_spec(idx_pk_max),'.','color',[0 0.8 0.1],'markersize',30)           
        plot(freqs(idx_high_trough), power_spec(idx_high_trough), '.','color',[0 0.4 1],'markersize',30)       
        plot(freqs(idx_low_trough), power_spec(idx_low_trough), '.','color',[0 0.4 1],'markersize',30)       
        text(10000, 0.8, strcat(num2str(round(pow)), ' uV^2'))

        Powers_Tecta(p)    = pow;
        Peak_freq_Tecta(p) = freqs(idx_pk_max); 

        ylim([0 1])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        if p == 1
            ylabel('Power Density $(\mu V^2 / Hz)$','Interpreter','latex','fontsize',14)
        end   
end

set(gcf,'Position',[0 100 1200 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_Tecta_Spectrograms_PowerSpecs.pdf'), 'ContentType', 'vector');

%% - Plot trace, spectrogram, power spectrum -> VGLUT3
clf
for p = 1:size(Data_VGLUT,2)
    subplot_tight(7,size(Data_VGLUT,2), p)
        signal      = Data_VGLUT(p).OCT{1,1};    
        idx_baseline_start  = (Data_VGLUT(p).laser_onset - 0.5)/0.01;
        idx_baseline_end    = (Data_VGLUT(p).laser_onset)/0.01;
        baseline    = mean(signal(idx_baseline_start : idx_baseline_end));
        moving_avg  = movmean(signal,200);        
        signal_z    = signal - moving_avg;
        signal_b    = signal - baseline;

        plot(Data_VGLUT(p).t_OCT, signal_z,'color','k','linewidth',2)
        text((Data_VGLUT(p).laser_onset - 0.85), 55, Data_VGLUT(p).filename,'interpreter','none','fontsize',8)
        title('VGLUT3')
        xlim([(Data_VGLUT(p).laser_onset - 1) (Data_VGLUT(p).laser_onset +3)])
        ylim([-70 70])

    subplot_tight(7, size(Data_VGLUT,2), [(p+3) (p+9)])
        window   = 70;      % 50
        noverlap = window-1;      % 75% of window length; 30
        nfft     = 10000;   % # of frequencies
        Fs       = 100000;  % Sampling rate

        spectrogram(signal_b ,window ,noverlap ,nfft ,Fs ,'yaxis'); % Short-time fourier transform
        colormap(flipud(turbo))
        clim([-45 0])
        ylim([0 15])
        xlim([(Data_VGLUT(p).laser_onset - 1) (Data_VGLUT(p).laser_onset +2)])

    subplot_tight(7, size(Data_VGLUT,2), [(p+12) (p+18)])
        t_start     = Data_VGLUT(p).laser_onset; % in ms
        t_end       = (Data_VGLUT(p).laser_onset) + 2;
        Fs          = 100000;
        [power_spec, freqs] = pspectrum(Data_VGLUT(p).OCT{1,1}(t_start/0.01 : t_end/0.01), Fs);
        plot(freqs, power_spec, 'k','linewidth',2)
        hold on

        [troughs, troughlocs]   = findpeaks(-power_spec); 
        idx_first_trough        = troughlocs(1);
        if p ==1
            idx_first_trough        = troughlocs(2);
        end
        [pk_max, idx_pk_max]    = findpeaks(power_spec(idx_first_trough:end), 'SortStr', 'descend','NPeaks',1);
        idx_pk_max              = idx_pk_max + idx_first_trough;
        [~ , idx_high_troughs]  = findpeaks(-power_spec(idx_pk_max:end));
        idx_high_trough         = idx_pk_max + idx_high_troughs(1);
        [~ , idx_low_troughs]   = findpeaks(-power_spec(1:idx_pk_max));
        idx_low_trough          = idx_low_troughs(end);
        pow = sum(power_spec(idx_low_trough:idx_high_trough)) * freqs(2);

        plot(freqs(1:idx_first_trough), power_spec(1:idx_first_trough), 'r','linewidth',2)    
        plot(freqs(idx_pk_max), power_spec(idx_pk_max),'.','color',[0 0.8 0.1],'markersize',30)           
        plot(freqs(idx_high_trough), power_spec(idx_high_trough), '.','color',[0 0.4 1],'markersize',30)       
        plot(freqs(idx_low_trough), power_spec(idx_low_trough), '.','color',[0 0.4 1],'markersize',30)       
        text(10000, 10, strcat(num2str(round(pow)), ' uV^2'),'fontsize',14)

        Powers_VGLUT3(p)    = pow;
        Peak_freq_VGLUT3(p) = freqs(idx_pk_max); 

        ylim([0 14])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        if p == 1
            ylabel('Power Density $(\mu V^2 / Hz)$','Interpreter','latex','fontsize',14)
        end           
        
end

set(gcf,'Position',[0 100 900 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_VGLUT3_Spectrograms_PowerSpecs.pdf'), 'ContentType', 'vector');

%% - Spectrogram - WT    - Band Integral
clf
% Input the frequency band over which to integrate:
integration_band = [5000 12000];

% Convert the frequency band into the corresponding array indices:
[~, idx_band_min] = (min(abs(freqs - integration_band(1))));
[~, idx_band_max] = (min(abs(freqs - integration_band(2))));

for p = 1:size(Data_WT,2)
    subplot_tight(7,size(Data_WT,2), p)
        signal      = Data_WT(p).OCT{1,1};    
        idx_baseline_start  = (Data_WT(p).laser_onset - 0.5)/0.01;
        idx_baseline_end    = (Data_WT(p).laser_onset)/0.01;
        baseline    = mean(signal(idx_baseline_start : idx_baseline_end));
        moving_avg  = movmean(signal,200);        
        signal_z    = signal - moving_avg;
        signal_b    = signal - baseline;

        plot(Data_WT(p).t_OCT, signal_z,'color','k','linewidth',2)
        text((Data_WT(p).laser_onset - 0.85), 55, Data_WT(p).filename,'interpreter','none','fontsize',8)
        title('WT')
        xlim([(Data_WT(p).laser_onset - 1) (Data_WT(p).laser_onset +3)])
        ylim([-70 70])

    subplot_tight(7, size(Data_WT,2), [(p+5) (p+15)])
        window   = 70;      % 50
        noverlap = window-1;      % 75% of window length; 30
        nfft     = 10000;   % # of frequencies
        Fs       = 100000;  % Sampling rate

        spectrogram(signal_b ,window ,noverlap ,nfft ,Fs ,'yaxis'); % Short-time fourier transform
        colormap(flipud(turbo))
        clim([-45 0])
        ylim([0 15])
        xlim([(Data_WT(p).laser_onset - 1) (Data_WT(p).laser_onset +2)])        

    subplot_tight(7, size(Data_WT,2), [(p+20) (p+30)])
        t_start     = Data_WT(p).laser_onset; % in ms
        t_end       = (Data_WT(p).laser_onset) + 2;
        Fs          = 100000;
        [power_spec, freqs] = pspectrum(Data_WT(p).OCT{1,1}(t_start/0.01 : t_end/0.01), Fs, 'power');
        plot(freqs, power_spec, 'k','linewidth',2)
        ylim([0 1])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        hold on

        %[troughs, troughlocs]   = findpeaks(-power_spec); 
        %idx_first_trough        = troughlocs(1);
        [pk_max, idx_pk_max]    = findpeaks(power_spec, 'SortStr', 'descend','NPeaks',1);
        %[~ , idx_high_troughs]  = findpeaks(-power_spec(idx_pk_max:end));
        %idx_high_trough         = idx_pk_max + idx_high_troughs(1);
        %[~ , idx_low_troughs]   = findpeaks(-power_spec(1:idx_pk_max));
        %idx_low_trough          = idx_low_troughs(end);

        pow = sum(power_spec(idx_band_min:idx_band_max)) * freqs(2);

        %plot(freqs(1:idx_first_trough), power_spec(1:idx_first_trough), 'r','linewidth',2)
        % Put a point at the peak:
        plot(freqs(idx_pk_max), power_spec(idx_pk_max),'.','color',[0 0.8 0.1],'markersize',30)           
        %plot(freqs(idx_band_min), power_spec(idx_band_min), '.','color',[0 0.4 1],'markersize',30)       
        %plot(freqs(idx_band_max), power_spec(idx_band_max), '.','color',[0 0.4 1],'markersize',30)       
        text(10000, 0.8, strcat(num2str(round(pow)), ' uV^2'),'fontsize',14)

        patch([freqs(idx_band_min); freqs(idx_band_min:idx_band_max); freqs(idx_band_max)],...
            [0; power_spec(idx_band_min:idx_band_max); 0],'r','facealpha',0.2,'linestyle','none')
        
        Powers_WT(p)    = pow;
        Peak_freq_WT(p) = freqs(idx_pk_max); 

        ylim([0 1])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        if p == 1
            ylabel('Power Density $(\mu V^2 / Hz)$','Interpreter','latex','fontsize',14)
        end    
        
end

set(gcf,'Position',[0 100 1500 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_WT_Spectrogram_BandIntegral.pdf'), 'ContentType', 'vector');

%% - Spectrogram - Tecta - Band Integral
clf
for p = 1:size(Data_Tecta,2)
    subplot_tight(7,size(Data_Tecta,2), p)
        signal      = Data_Tecta(p).OCT{1,1};    
        idx_baseline_start  = (Data_Tecta(p).laser_onset - 0.5)/0.01;
        idx_baseline_end    = (Data_Tecta(p).laser_onset)/0.01;
        baseline    = mean(signal(idx_baseline_start : idx_baseline_end));
        moving_avg  = movmean(signal,200);        
        signal_z    = signal - moving_avg;
        signal_b    = signal - baseline;

        plot(Data_Tecta(p).t_OCT, signal_z,'color','k','linewidth',2)
        text((Data_Tecta(p).laser_onset - 0.85), 55, Data_Tecta(p).filename,'interpreter','none','fontsize',8)
        title('Tecta')
        xlim([(Data_Tecta(p).laser_onset - 1) (Data_Tecta(p).laser_onset +3)])
        ylim([-70 70])

    subplot_tight(7, size(Data_Tecta,2), [(p+4) (p+12)])
        window   = 70;      % 50
        noverlap = window-1;      % 75% of window length; 30
        nfft     = 10000;   % # of frequencies
        Fs       = 100000;  % Sampling rate

        spectrogram(signal_b ,window ,noverlap ,nfft ,Fs ,'yaxis'); % Short-time fourier transform
        colormap(flipud(turbo))
        clim([-45 0])
        ylim([0 15])
        xlim([(Data_Tecta(p).laser_onset - 1) (Data_Tecta(p).laser_onset +2)])

    subplot_tight(7, size(Data_Tecta,2), [(p+16) (p+24)])
        t_start     = Data_Tecta(p).laser_onset; % in ms
        t_end       = (Data_Tecta(p).laser_onset) + 2;
        Fs          = 100000;
        [power_spec, freqs] = pspectrum(Data_Tecta(p).OCT{1,1}(t_start/0.01 : t_end/0.01), Fs);
        plot(freqs, power_spec, 'k','linewidth',2)
        hold on

        %[troughs, troughlocs]   = findpeaks(-power_spec); 
        %idx_first_trough        = troughlocs(1);
        [pk_max, idx_pk_max]    = findpeaks(power_spec, 'SortStr', 'descend','NPeaks',1);
        %[~ , idx_high_troughs]  = findpeaks(-power_spec(idx_pk_max:end));
        %idx_high_trough         = idx_pk_max + idx_high_troughs(1);
        %[~ , idx_low_troughs]   = findpeaks(-power_spec(1:idx_pk_max));
        %idx_low_trough          = idx_low_troughs(end);

        pow = sum(power_spec(idx_band_min:idx_band_max)) * freqs(2);

        %plot(freqs(1:idx_first_trough), power_spec(1:idx_first_trough), 'r','linewidth',2)
        % Put a point at the peak:
        plot(freqs(idx_pk_max), power_spec(idx_pk_max),'.','color',[0 0.8 0.1],'markersize',30)           
        %plot(freqs(idx_band_min), power_spec(idx_band_min), '.','color',[0 0.4 1],'markersize',30)       
        %plot(freqs(idx_band_max), power_spec(idx_band_max), '.','color',[0 0.4 1],'markersize',30)       
        text(10000, 0.8, strcat(num2str(round(pow)), ' uV^2'),'fontsize',14)

        patch([freqs(idx_band_min); freqs(idx_band_min:idx_band_max); freqs(idx_band_max)],...
            [0; power_spec(idx_band_min:idx_band_max); 0],'r','facealpha',0.2,'linestyle','none')
        
        Powers_Tecta(p)    = pow;
        Peak_freq_Tecta(p) = freqs(idx_pk_max); 

        ylim([0 1])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        if p == 1
            ylabel('Power Density $(\mu V^2 / Hz)$','Interpreter','latex','fontsize',14)
        end   
end

set(gcf,'Position',[0 100 1200 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_Tecta_Spectrograms_BandIntegral.pdf'), 'ContentType', 'vector');

%% - Spectrogram - VGLUT - Band Integral
clf
for p = 1:size(Data_VGLUT,2)
    subplot_tight(7,size(Data_VGLUT,2), p)
        signal      = Data_VGLUT(p).OCT{1,1};    
        idx_baseline_start  = (Data_VGLUT(p).laser_onset - 0.5)/0.01;
        idx_baseline_end    = (Data_VGLUT(p).laser_onset)/0.01;
        baseline    = mean(signal(idx_baseline_start : idx_baseline_end));
        moving_avg  = movmean(signal,200);        
        signal_z    = signal - moving_avg;
        signal_b    = signal - baseline;

        plot(Data_VGLUT(p).t_OCT, signal_z,'color','k','linewidth',2)
        text((Data_VGLUT(p).laser_onset - 0.85), 55, Data_VGLUT(p).filename,'interpreter','none','fontsize',8)
        title('VGLUT3')
        xlim([(Data_VGLUT(p).laser_onset - 1) (Data_VGLUT(p).laser_onset +3)])
        ylim([-70 70])

    subplot_tight(7, size(Data_VGLUT,2), [(p+3) (p+9)])
        window   = 70;      % 50
        noverlap = window-1;      % 75% of window length; 30
        nfft     = 10000;   % # of frequencies
        Fs       = 100000;  % Sampling rate

        spectrogram(signal_b ,window ,noverlap ,nfft ,Fs ,'yaxis'); % Short-time fourier transform
        colormap(flipud(turbo))
        clim([-45 0])
        ylim([0 15])
        xlim([(Data_VGLUT(p).laser_onset - 1) (Data_VGLUT(p).laser_onset +2)])

    subplot_tight(7, size(Data_VGLUT,2), [(p+12) (p+18)])
        t_start     = Data_VGLUT(p).laser_onset; % in ms
        t_end       = (Data_VGLUT(p).laser_onset) + 2;
        Fs          = 100000;
        [power_spec, freqs] = pspectrum(Data_VGLUT(p).OCT{1,1}(t_start/0.01 : t_end/0.01), Fs);
        plot(freqs, power_spec, 'k','linewidth',2)
        hold on

        %[troughs, troughlocs]   = findpeaks(-power_spec); 
        %idx_first_trough        = troughlocs(1);
        
        if p == 1
            [pk_max, idx_pk_max]    = findpeaks(power_spec, 'SortStr', 'descend','NPeaks',2);
            pk_max = pk_max(2);
            idx_pk_max = idx_pk_max(2);
        else
            [pk_max, idx_pk_max]    = findpeaks(power_spec, 'SortStr', 'descend','NPeaks',1);
        end
        %[~ , idx_high_troughs]  = findpeaks(-power_spec(idx_pk_max:end));
        %idx_high_trough         = idx_pk_max + idx_high_troughs(1);
        %[~ , idx_low_troughs]   = findpeaks(-power_spec(1:idx_pk_max));
        %idx_low_trough          = idx_low_troughs(end);

        pow = sum(power_spec(idx_band_min:idx_band_max)) * freqs(2);

        %plot(freqs(1:idx_first_trough), power_spec(1:idx_first_trough), 'r','linewidth',2)
        % Put a point at the peak:
        plot(freqs(idx_pk_max), power_spec(idx_pk_max),'.','color',[0 0.8 0.1],'markersize',30)           
        %plot(freqs(idx_band_min), power_spec(idx_band_min), '.','color',[0 0.4 1],'markersize',30)       
        %plot(freqs(idx_band_max), power_spec(idx_band_max), '.','color',[0 0.4 1],'markersize',30)       
        text(12000, 8, strcat(num2str(round(pow)), ' uV^2'),'fontsize',14)

        patch([freqs(idx_band_min); freqs(idx_band_min:idx_band_max); freqs(idx_band_max)],...
            [0; power_spec(idx_band_min:idx_band_max); 0],'r','facealpha',0.2,'linestyle','none')
        
        Powers_VGLUT3(p)    = pow;
        Peak_freq_VGLUT3(p) = freqs(idx_pk_max); 

        ylim([0 14])
        xlim([0 15000])
        xlabel('Freq (Hz)')
        if p == 1
            ylabel('Power Density $(\mu V^2 / Hz)$','Interpreter','latex','fontsize',14)
        end   
end

set(gcf,'Position',[0 100 900 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_VGLUT3_Spectrograms_BandIntegral.pdf'), 'ContentType', 'vector');

%% - Histogram of powers for each mutant
clf
plot(0.2+0.4.*(rand(1,length(Powers_WT))),log10(Powers_WT),'.k','markersize',40)
hold on
plot(1.0+0.4.*(rand(1,length(Powers_Tecta))),log10(Powers_Tecta),'.k','markersize',40)
plot(1.9+0.4.*(rand(1,length(Powers_VGLUT3))),log10(Powers_VGLUT3),'.k','markersize',40)

% [H, p] = unpaired_ttest(Powers_WT, Powers_Tecta, 1, 4, [0.2 1]);
% [H, p] = unpaired_ttest(Powers_WT, Powers_VGLUT3, 1, 4.8, [0.2 1.9]);
% [H, p] = unpaired_ttest(Powers_Tecta, Powers_VGLUT3, 1, 4.5, [1 1.9]);

[p] = KruskalWallisPlot(Powers_WT', Powers_Tecta', 1, 4, [0.2 1]);
[p] = KruskalWallisPlot(Powers_WT', Powers_VGLUT3', 1, 4.8, [0.2 1.9]);
[p] = KruskalWallisPlot(Powers_Tecta', Powers_VGLUT3', 1, 4.5, [1 1.9]);

xlim([0 2.4])
xticks([0.2 1.1 2])
set(gca,'fontsize',18,'xticklabel',{'WT', 'Tecta', 'VGLUT3'})
ylabel('Log (Power uV^2)')
ylim([0 5])
set(gcf,'Position',[800 600 500 400])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_C_Mutants_OCT_RingingPower_Histogram.pdf'), 'ContentType', 'vector');

%% - Make OCT tuning curves - WT

folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";

clf
for p = 1:size(Data_WT,2)

    [Freqs, Ampls, Data, Gain ] = OCTtuningCurve(folder, strcat('0',num2str(Data_WT(p).Date)), Data_WT(p).Mouse, Data_WT(p).file_tuning);
    
subplot_tight(2,size(Data_WT,2),p)
    plot(Freqs{:,1}, log10(Data{:,:}),'linewidth',2)
    text(0, 2.4, Data_WT(p).file_tuning,'interpreter','none','fontsize',8)
    ylim([-1.5 2.5])
    title('WT')
    set(gca,'ColorOrderIndex',1,'fontsize',10)
    legend('0 dB','10 dB','20 dB','30 dB','40 dB','50 dB','60 dB','70 dB','Location','southwest')
    ylabel('log (displacement)')

subplot_tight(2,size(Data_WT,2),p+size(Data_WT,2))    
    plot(Freqs{:,1}, (Gain(:,1)),'linewidth',2,'color',[0.9 0.9 0.9])
    hold on
    plot(Freqs{:,1}, (Gain(:,2:end)),'linewidth',2)
    legend('0 dB','10 dB','20 dB','30 dB','40 dB','50 dB','60 dB','70 dB','Location','northeast')
    set(gca, 'fontsize',10)
    ylim([0 3.5e4])
    ylabel('Gain (displ / Pressure in uPa)')
    %ylabel('Gain (Displ / Pressure in uPa)')
end

%TuningPeaks_WT = [NaN; 9000; NaN; 5000; 6000];
TuningPeaks_WT = [9000; 9000; 6000; 5000; 6000];

%% - Make OCT tuning curves - Tecta

folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";

clf
for p = 1:size(Data_Tecta,2)

    [Freqs, Ampls, Data, Gain ] = OCTtuningCurve(folder, strcat('0',num2str(Data_Tecta(p).Date)), Data_Tecta(p).Mouse, Data_Tecta(p).file_tuning);
    
subplot_tight(2,size(Data_Tecta,2),p)
    plot(Freqs{:,1}, log10(Data{:,:}),'linewidth',2)
    text(0, 2.4, Data_Tecta(p).file_tuning,'interpreter','none','fontsize',8)
    ylim([-1.5 2.5])
    title('Tecta')
    set(gca,'ColorOrderIndex',1,'fontsize',10)
    legend('0 dB','10 dB','20 dB','30 dB','40 dB','50 dB','60 dB','70 dB','Location','southwest')
    ylabel('log (displacement)')

subplot_tight(2,size(Data_Tecta,2),p+size(Data_Tecta,2))    
    plot(Freqs{:,1}, (Gain(:,1)),'linewidth',2,'color',[0.9 0.9 0.9])
    hold on
    plot(Freqs{:,1}, (Gain(:,2:end)),'linewidth',2)
    legend('0 dB','10 dB','20 dB','30 dB','40 dB','50 dB','60 dB','70 dB','Location','northeast')
    set(gca, 'fontsize',10)
    ylim([0 3e3])
    ylabel('Gain (displ / Pressure in uPa)')
    %ylabel('Gain (Displ / Pressure in uPa)')
end

TuningPeaks_Tecta = [4000; 3500; 4000; 4000];

%% - Make OCT tuning curves - VGLUT3

folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Oghalai/";

clf
for p = 1:size(Data_VGLUT,2)

    [Freqs, Ampls, Data, Gain ] = OCTtuningCurve(folder, strcat('0',num2str(Data_VGLUT(p).Date)), Data_VGLUT(p).Mouse, Data_VGLUT(p).file_tuning);
    
subplot_tight(2,size(Data_VGLUT,2),p)
    plot(Freqs{:,1}, log10(Data{:,:}),'linewidth',2)
    text(0, 2.4, Data_VGLUT(p).file_tuning,'interpreter','none','fontsize',8)
    ylim([-1.5 2.5])
    title('VGLUT3')
    set(gca,'ColorOrderIndex',1,'fontsize',10)
    legend('0 dB','10 dB','20 dB','30 dB','40 dB','50 dB','60 dB','70 dB','Location','southwest')
    ylabel('log (displacement)')

subplot_tight(2,size(Data_VGLUT,2),p+size(Data_VGLUT,2))    
    plot(Freqs{:,1}, (Gain(:,1)),'linewidth',2,'color',[0.9 0.9 0.9])
    hold on
    plot(Freqs{:,1}, (Gain(:,2:end)),'linewidth',2)
    legend('0 dB','10 dB','20 dB','30 dB','40 dB','50 dB','60 dB','70 dB','Location','northeast')
    set(gca, 'fontsize',10)
    ylim([0 45e4])
    ylabel('Gain (displ / Pressure in uPa)')
    %ylabel('Gain (Displ / Pressure in uPa)')
end

TuningPeaks_VGLUT3 = [8000; 9500; 10000];

%% - Plot OCT Sound Tuning Freq vs Laser PSD peak freq
clf

plot(TuningPeaks_WT, Peak_freq_WT,'.k','markersize',40)
    hold on
plot(TuningPeaks_Tecta, Peak_freq_Tecta,'.b','markersize',40)
plot(TuningPeaks_VGLUT3, Peak_freq_VGLUT3,'.r','markersize',40)
    text(2100, 9800,'WT','fontsize',16)
    text(2100, 9400,'Tecta','fontsize',16,'color','b')
    text(2100, 9000,'VGLUT3','fontsize',16,'color','r')
    xlabel('Peak on OCT tuning curve (Hz)')
    ylabel('Frequency of peak on PSD of laser response (Hz)')
    ylim([2000 10000])
    xlim([2000 11000])    
    
    set(gca,'fontsize',16)
    set(gcf,'Position',[800 500 600 500])

% - Calculate the linear regression
AllTuningPeaks  = [TuningPeaks_WT; TuningPeaks_Tecta; TuningPeaks_VGLUT3];
AllLaserPeaks   = [Peak_freq_WT'; Peak_freq_Tecta'; Peak_freq_VGLUT3'];

pfit = polyfit(AllTuningPeaks, AllLaserPeaks, 1);

x       = linspace(3000, 10000);
yfit    = x.*pfit(1) + pfit(2);

plot(x, yfit,'color',[0.8 0.8 0.8],'linewidth',2)
text(8500,5000,strcat('slope = ', num2str(pfit(1))),'fontsize',16)
text(8500,4500,strcat('y-int = ', num2str(pfit(2))),'fontsize',16)

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_B_Mutants_TuningF_Vs_LaserPSD_peak.pdf'), 'ContentType', 'vector');

%% - Calculate and plot oCAP1 and oCAP2 for WT, Tecta, VGLUT3
close all
% - Wild Type - %
for j = 1: size(Data_WT,2)
    subplot(size(Data_WT,2),1,j)
        signal  = Data_WT(j).CAPs(1,:);
        t       = Data_WT(j).t_CAP;
        plot(t, signal,'k','linewidth',1.5)
        hold on

        xlimvec = [(Data_WT(j).laser_onset)-1 (Data_WT(j).laser_onset)+6];
        xlim(xlimvec)

        oCAP1_min_start = (Data_WT(j).laser_onset)-0.1;
        oCAP1_min_end   = (Data_WT(j).laser_onset)+0.2;
        oCAP1_max_start = (Data_WT(j).laser_onset)+0.1;
        oCAP1_max_end   = (Data_WT(j).laser_onset)+0.4;

        oCAP2_min_start = (Data_WT(j).laser_onset)+0.6;
        oCAP2_min_end   = (Data_WT(j).laser_onset)+1.6;
        oCAP2_max_start = (Data_WT(j).laser_onset)+1.5;
        oCAP2_max_end   = (Data_WT(j).laser_onset)+2.1;

        % Convert the times in ms to indices of the vector t
        [~, oCAP1_min_start_idx]  = findNearest(Data_WT(j).t_CAP, oCAP1_min_start);
        [~, oCAP1_min_end_idx]    = findNearest(Data_WT(j).t_CAP, oCAP1_min_end);
        [~, oCAP1_max_start_idx]  = findNearest(Data_WT(j).t_CAP, oCAP1_max_start);
        [~, oCAP1_max_end_idx]    = findNearest(Data_WT(j).t_CAP, oCAP1_max_end);
        
        % Convert the times in ms to indices of the vector t
        [~, oCAP2_min_start_idx]  = findNearest(Data_WT(j).t_CAP, oCAP2_min_start);
        [~, oCAP2_min_end_idx]    = findNearest(Data_WT(j).t_CAP, oCAP2_min_end);
        [~, oCAP2_max_start_idx]  = findNearest(Data_WT(j).t_CAP, oCAP2_max_start);
        [~, oCAP2_max_end_idx]    = findNearest(Data_WT(j).t_CAP, oCAP2_max_end);
        
         % Find peak 1
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
        WT_oCAP1_Ampl(j) = Ampl;
        % Markers for peak 1
        plot(t(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t(idx_max),signal(idx_max),'.g','markersize',20)
        % Find peak 2
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
        WT_oCAP2_Ampl(j) = Ampl;
        % Markers for peak 2
        plot(t(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t(idx_max),signal(idx_max),'.g','markersize',20)
        text(Data_WT(j).laser_onset-0.5, 40, Data_WT(j).filename,'interpreter','none')

        if j == 1
            title('Wild Type')
        end

        ylim([-60 60])
end

set(gcf,'Position',[1400 0 400 1000])

% % % % % % % % % % % % % % % % % % % % 
% - Tecta - %
figure
for j = 1: size(Data_Tecta,2)
    subplot(size(Data_Tecta,2),1,j)
        signal  = Data_Tecta(j).CAPs(1,:);
        t       = Data_Tecta(j).t_CAP;
        plot(t, signal,'k','linewidth',1.5)
        hold on

        xlimvec = [(Data_Tecta(j).laser_onset)-1 (Data_Tecta(j).laser_onset)+6];
        xlim(xlimvec)

        oCAP1_min_start = (Data_Tecta(j).laser_onset)-0.1;
        oCAP1_min_end   = (Data_Tecta(j).laser_onset)+0.2;
        oCAP1_max_start = (Data_Tecta(j).laser_onset)+0.1;
        oCAP1_max_end   = (Data_Tecta(j).laser_onset)+0.4;

        oCAP2_min_start = (Data_Tecta(j).laser_onset)+0.6;
        oCAP2_min_end   = (Data_Tecta(j).laser_onset)+1.6;
        oCAP2_max_start = (Data_Tecta(j).laser_onset)+1.5;
        oCAP2_max_end   = (Data_Tecta(j).laser_onset)+2.1;

        % Convert the times in ms to indices of the vector t
        [~, oCAP1_min_start_idx]  = findNearest(Data_Tecta(j).t_CAP, oCAP1_min_start);
        [~, oCAP1_min_end_idx]    = findNearest(Data_Tecta(j).t_CAP, oCAP1_min_end);
        [~, oCAP1_max_start_idx]  = findNearest(Data_Tecta(j).t_CAP, oCAP1_max_start);
        [~, oCAP1_max_end_idx]    = findNearest(Data_Tecta(j).t_CAP, oCAP1_max_end);
        
        % Convert the times in ms to indices of the vector t
        [~, oCAP2_min_start_idx]  = findNearest(Data_Tecta(j).t_CAP, oCAP2_min_start);
        [~, oCAP2_min_end_idx]    = findNearest(Data_Tecta(j).t_CAP, oCAP2_min_end);
        [~, oCAP2_max_start_idx]  = findNearest(Data_Tecta(j).t_CAP, oCAP2_max_start);
        [~, oCAP2_max_end_idx]    = findNearest(Data_Tecta(j).t_CAP, oCAP2_max_end);
        
         % Find peak 1
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
        Tecta_oCAP1_Ampl(j) = Ampl;
        % Markers for peak 1
        plot(t(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t(idx_max),signal(idx_max),'.g','markersize',20)
        % Find peak 2
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
        Tecta_oCAP2_Ampl(j) = Ampl;
        % Markers for peak 2
        plot(t(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t(idx_max),signal(idx_max),'.g','markersize',20)
        text(Data_Tecta(j).laser_onset-0.5, 40, Data_Tecta(j).filename,'interpreter','none')

        if j == 1
            title('Tecta')
        end

        ylim([-60 60])
end

set(gcf,'Position',[900 0 400 1000])


% % % % % % % % % % % % % % % % % % % % 
% - VGLUT3 - %
figure
for j = 1: size(Data_VGLUT,2)
    subplot(size(Data_VGLUT,2),1,j)
        signal  = Data_VGLUT(j).CAPs(1,:);
        t       = Data_VGLUT(j).t_CAP;
        plot(t, signal,'k','linewidth',1.5)
        hold on

        xlimvec = [(Data_VGLUT(j).laser_onset)-1 (Data_VGLUT(j).laser_onset)+6];
        xlim(xlimvec)

        oCAP1_min_start = (Data_VGLUT(j).laser_onset)-0.1;
        oCAP1_min_end   = (Data_VGLUT(j).laser_onset)+0.25;
        oCAP1_max_start = (Data_VGLUT(j).laser_onset)+0.1;
        oCAP1_max_end   = (Data_VGLUT(j).laser_onset)+0.5;

        oCAP2_min_start = (Data_VGLUT(j).laser_onset)+0.6;
        oCAP2_min_end   = (Data_VGLUT(j).laser_onset)+1.6;
        oCAP2_max_start = (Data_VGLUT(j).laser_onset)+1.5;
        oCAP2_max_end   = (Data_VGLUT(j).laser_onset)+2.1;

        % Convert the times in ms to indices of the vector t
        [~, oCAP1_min_start_idx]  = findNearest(Data_VGLUT(j).t_CAP, oCAP1_min_start);
        [~, oCAP1_min_end_idx]    = findNearest(Data_VGLUT(j).t_CAP, oCAP1_min_end);
        [~, oCAP1_max_start_idx]  = findNearest(Data_VGLUT(j).t_CAP, oCAP1_max_start);
        [~, oCAP1_max_end_idx]    = findNearest(Data_VGLUT(j).t_CAP, oCAP1_max_end);
        
        % Convert the times in ms to indices of the vector t
        [~, oCAP2_min_start_idx]  = findNearest(Data_VGLUT(j).t_CAP, oCAP2_min_start);
        [~, oCAP2_min_end_idx]    = findNearest(Data_VGLUT(j).t_CAP, oCAP2_min_end);
        [~, oCAP2_max_start_idx]  = findNearest(Data_VGLUT(j).t_CAP, oCAP2_max_start);
        [~, oCAP2_max_end_idx]    = findNearest(Data_VGLUT(j).t_CAP, oCAP2_max_end);
        
         % Find peak 1
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
        VGLUT3_oCAP1_Ampl(j) = Ampl;
        % Markers for peak 1
        plot(t(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t(idx_max),signal(idx_max),'.g','markersize',20)
        % Find peak 2
        [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
        VGLUT3_oCAP2_Ampl(j) = Ampl;
        % Markers for peak 2
        plot(t(idx_min),signal(idx_min),'.r','markersize',20)
        plot(t(idx_max),signal(idx_max),'.g','markersize',20)
        text(Data_VGLUT(j).laser_onset-0.5, 40, Data_VGLUT(j).filename,'interpreter','none')

        if j == 1
            title('VGLUT3')
        end

        ylim([-60 60])
end

set(gcf,'Position',[500 0 400 1000])

%% - Plot CAP summary data
close all
clf
subplot(1,2,1)
plot(0.4.*(rand(1,length(WT_oCAP1_Ampl))),WT_oCAP1_Ampl,'.k','markersize',40)
hold on
plot(0.4+0.35.*(rand(1,length(Tecta_oCAP1_Ampl))),Tecta_oCAP1_Ampl,'.r','markersize',40)
xlim([0 0.8])
xticks([0.2 0.6])
title('oCAP1')
ylabel('Peak-peak (uV)')
set(gca,'fontsize',16,'xticklabel',{'WT', 'Tecta'})
%[H, p] = unpaired_ttest(WT_oCAP1_Ampl, Tecta_oCAP1_Ampl, 1, 100, [0.2 0.6]);
[p] = KruskalWallisPlot(WT_oCAP1_Ampl', Tecta_oCAP1_Ampl', 1, 100, [0.2 0.6]);
ylim([0 120])

subplot(1,2,2)
plot(0.4.*(rand(1,length(WT_oCAP2_Ampl))),WT_oCAP2_Ampl,'.k','markersize',40)
hold on
plot(0.4+0.35.*(rand(1,length(Tecta_oCAP2_Ampl))),Tecta_oCAP2_Ampl,'.r','markersize',40)
xlim([0 0.8])
xticks([0.2 0.6])
title('oCAP2')
ylabel('Peak-peak (uV)')
set(gca,'fontsize',16,'xticklabel',{'WT', 'Tecta'})
%[H, p] = unpaired_ttest(WT_oCAP2_Ampl, Tecta_oCAP2_Ampl, 1, 100, [0.2 0.6]);
[p] = KruskalWallisPlot(WT_oCAP2_Ampl', Tecta_oCAP2_Ampl', 1, 100, [0.2 0.6]);
ylim([0 120])

set(gcf,'Position',[1200 700 550 260])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, strcat('Figures_OCT_Fig2_C_Summary_oCAPs_WT_Tecta.pdf'), 'ContentType', 'vector');




