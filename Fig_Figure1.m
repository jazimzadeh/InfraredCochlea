%% - Part B - ABR and CAP in response to click
% Use 2023_1_26
clf
clear
exptdate = '2023_1_26';
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
path = strcat(folder,exptdate,'/');

cd(path)
clf

RunOne          = readmatrix(strcat('Run1.csv'));

Fs          = 25000;                    % Sampling rate in Hz
tstep       = 1/Fs;
t           = [0:(length(RunOne)-49)]; % # of time points
t           = t.*tstep;                 % Time vector
laser_dur   = 0.150;                    % Laser pulse duration in ms

subplot(1,2,1)
offset = 0;
for j = 2:2:19
   plot(t.*(1e3),(RunOne(j,49:end) + offset).*0.001,'linewidth',1)
   hold on
   offset = offset - (3E3);
end
grid on
ylim([-25 5])
xlim([0 10])
title(strcat(exptdate, '- Run 1 ABR'), 'Interpreter', 'none')
xlabel('Time (ms)')
ylabel('Voltage uV')
set(gca, 'fontsize',13)

%legend('90 dB', '80 dB', '70 dB', '60 dB', '50 dB','40 dB','30 dB', '20 dB', '10 dB','location','southeast')

subplot(1,2,2)
offset = 0;
for j = 3:2:19
   plot(t.*(1e3),(RunOne(j,49:end) + offset).*0.001,'linewidth',1)
   hold on
   offset = offset - (30E3);
end
grid on
ylim([-250 50])
xlim([0 10])
title(strcat(exptdate, '- Run1 CAP'), 'Interpreter', 'none')
xlabel('Time (ms)')
ylabel('Voltage uV')
set(gca, 'fontsize',13)

legend('90 dB', '80 dB', '70 dB', '60 dB', '50 dB','40 dB','30 dB', '20 dB', '10 dB','location','southeast')
set(gcf,'Position',[0 700 500 400])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figure1_PartB_Click_ABR_CAP.pdf', 'ContentType', 'vector');

%% - Part C - Break down CAP phases
% Again using the data from 1/26/23, as imported above
cd(path)
clf

Run         = readmatrix(strcat('Run1-0-1-2-1.csv'));
RunAll      = mean(Run,1).*1000000;
RunEven     = mean(Run(1:2:end,:).*1000000,1);
RunOdd      = mean(Run(2:2:end,:).*1000000,1);

t = linspace(0,0.0107,2144).*1000;

plot(t,RunAll,'linewidth',1,'color','k')
hold on
plot(t,RunEven-60,'linewidth',1.,'color','b')
plot(t,RunOdd-60,'linewidth',1,'color','r')
plot(t,RunEven - RunAll-130,'linewidth',1,'color','b')
plot(t,RunOdd - RunAll-130,'linewidth',1,'color','r')


text(5,10,'Click Averaged','fontsize',14)
text(5,-50,'Each phase','fontsize',14)
text(3.5,-120,'Subtract average (Microphonic)','fontsize',14)


title('CAPs - Run 1 (Click) w/microphonic, Run 7 (Capella)')
xlabel('Time (ms)')
ylabel('\muV')
xlim([0 8])

set(gcf,'Position',[0 700 500 400])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figure1_PartC_CAP_phases_microphonic.pdf', 'ContentType', 'vector');

%% - Part E - Laser Response - 400 um fiber - oCAP1 and oCAP 2 vs power
% - 2023_3_19
clf
clear
exptdate = '2023_3_19';
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
path = strcat(folder,exptdate,'/');

% - Make 400 um oCAP1 and oCAP2 plots
Runs            = [12 13 14 15 16 17 18 19 20 21 22];
PowerPercent    = [100 80 60 50 40 30 20 15 17 14 14.5];

Target_Distance = 300;
Pulse_Length    = 100;
Fiber_Diameter  = 400;

cd(path)
% Make 3-D array containing all the CSV data
for i = 1:length(Runs)
    s(i,:,:) = readmatrix(strcat('Run', num2str(Runs(i)), '.csv'));
end

t  = linspace(0,0.01098,269).*1000;

close all

% Settings for plot offset and y-limit 
offset_oCAP  = 120;
xlimvec      = [0 10];

% Define the times within which to find oCAP peak 1 amplitude
oCAP1_min_start   = 0.9; 
oCAP1_min_end     = 1.1;
oCAP1_max_start   = 1.1;
oCAP1_max_end     = 1.4;

% Define the times within which to find oCAP peak 2 amplitude
oCAP2_min_start   = 1.5; 
oCAP2_min_end     = 2.5;
oCAP2_max_start   = 2.0; 
oCAP2_max_end     = 3;

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

offset = 0;
offset_oCAP = 35;
for p = 1:size(Runs,2)    
    signal = squeeze(s(p,3,49:end)).*(1e-3) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on  
     % Find peak 1
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
    if isempty(Ampl)
       oCAP1_Ampl(p) = 0;
    else
        oCAP1_Ampl(p) = Ampl;
    end
    % Markers for peak 1
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Find peak 2
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
    if isempty(Ampl)
        oCAP2_Ampl(p) = 0;
    elseif Ampl < 0
        oCAP2_Ampl(p) = 0;
    else
        oCAP2_Ampl(p) = Ampl;
    end
    % Markers for peak 2
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    
    offset = offset + offset_oCAP;
    grid on   
    xlabel('ms')
    %ylim(ylimvec_aCAP)
    xlim(xlimvec)

    [Watts(p), ~, ~, EnergyDensity(p)] = CapellaEnergy(Fiber_Diameter, PowerPercent(p), Target_Distance, Pulse_Length);
end

set(gcf,'Position',[0 0 400 1000])



figure
plot(EnergyDensity, oCAP2_Ampl, '.k','markersize', 30)
hold on
plot(EnergyDensity, oCAP1_Ampl, '.r','markersize', 30)
xlabel('Radiant Energy (mJ/cm^2)')
ylabel('uV')
set(gca, 'fontsize',15)
set(gcf, 'Position',[600 500 700 400])
legend('oCAP2', 'oCAP1','location', 'northwest')
title(strcat(exptdate,' - Mouse 1 - 400 um Fiber'),'interpreter','none')
set(gcf,'Position',[400 800 300 300])
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figure1_PartE_oCAP_PowerCurve.pdf', 'ContentType', 'vector');

figure
plot(EnergyDensity, oCAP2_Ampl, '.k','markersize', 30)
hold on
plot(EnergyDensity, oCAP1_Ampl, '.r','markersize', 30)
xlabel('Radiant Energy (mJ/cm^2)')
ylabel('uV')
set(gca, 'fontsize',15)
set(gcf, 'Position',[600 500 700 400])
legend('oCAP2','oCAP1','location', 'northwest')
title(strcat(exptdate,' - Mouse 1 - 400 um Fiber'),'interpreter','none')
xlim([0 20])
set(gcf,'Position',[400 400 300 300])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figure1_PartE_oCAP_PowerCurve_inset.pdf', 'ContentType', 'vector');

%% - Supplement: 200 vs 400 um fiber - Plot 200 um fiber oCAP1 and oCAP2 vs power
% - 2023_3_19
close all

exptdate = '2023_3_19';
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
path = strcat(folder,exptdate,'/');

% - Make 400 um oCAP1 and oCAP2 plots
%Runs            = [12 13 14 15 16 17 18 19 20 21 22];
%PowerPercent    = [100 80 60 50 40 30 20 15 17 14 14.5];

Runs_200        = [2 3 4 5 6 7 8 9 10 ];
PowerPercent200 = [100 80 60 50 40 30 20 15 17 ];

Target_Distance = 300;
Pulse_Length    = 100;
Fiber_Diameter  = 200;

cd(path)
% Make 3-D array containing all the CSV data
for i = 1:length(Runs_200)
    s(i,:,:) = readmatrix(strcat('Run', num2str(Runs_200(i)), '.csv'));
end

t  = linspace(0,0.01098,269).*1000;

close all

% Settings for plot offset and y-limit 
offset_oCAP  = 120;
xlimvec      = [0 10];

% Define the times within which to find oCAP peak 1 amplitude
oCAP1_min_start   = 0.9; 
oCAP1_min_end     = 1.1;
oCAP1_max_start   = 1.1;
oCAP1_max_end     = 1.4;

% Define the times within which to find oCAP peak 2 amplitude
oCAP2_min_start   = 1.5; 
oCAP2_min_end     = 2.5;
oCAP2_max_start   = 2.0; 
oCAP2_max_end     = 3;

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

offset = 0;
offset_oCAP = 35;
for p = 1:size(Runs_200,2)    
    signal = squeeze(s(p,3,49:end)).*(1e-3) + offset;
    plot(t, signal, 'k','linewidth',2)
    hold on  
     % Find peak 1
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP1_min_start_idx, oCAP1_min_end_idx, oCAP1_max_start_idx, oCAP1_max_end_idx);
    if isempty(Ampl)
        oCAP1_Ampl200(p) = 0;
    else
        oCAP1_Ampl200(p) = Ampl;
    end
    % Markers for peak 1
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    % Find peak 2
    [Ampl, idx_min, idx_max] = CAP_Ampl(signal, oCAP2_min_start_idx, oCAP2_min_end_idx, oCAP2_max_start_idx, oCAP2_max_end_idx);
    if isempty(Ampl)
        oCAP2_Ampl200(p) = 0;
    elseif Ampl < 0
        oCAP2_Ampl200(p) = 0;
    else
        oCAP2_Ampl200(p) = Ampl;
    end
    % Markers for peak 2
    plot(t(idx_min),signal(idx_min),'.r','markersize',20)
    plot(t(idx_max),signal(idx_max),'.g','markersize',20)
    
    offset = offset + offset_oCAP;
    grid on   
    xlabel('ms')
    %ylim(ylimvec_aCAP)
    xlim(xlimvec)

    [Watts200(p), ~, ~, EnergyDensity200(p)] = CapellaEnergy(Fiber_Diameter, PowerPercent200(p), Target_Distance, Pulse_Length);
end

set(gcf,'Position',[0 0 400 1000])



figure
plot(EnergyDensity200, oCAP2_Ampl200, '.k','markersize', 30)
hold on
plot(EnergyDensity200, oCAP1_Ampl200, '.r','markersize', 30)
xlabel('Radiant Energy (mJ/cm^2)')
ylabel('uV')
set(gca, 'fontsize',15)
set(gcf, 'Position',[600 500 700 400])
legend('oCAP2', 'oCAP1','location', 'northwest')
title(strcat(exptdate,' - Mouse 1 - 200 um Fiber'),'interpreter','none')
set(gcf,'Position',[400 800 300 300])
%cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
%exportgraphics(gcf, 'Figure1_PartE_oCAP_PowerCurve.pdf', 'ContentType', 'vector');

figure
plot(EnergyDensity200, oCAP2_Ampl200, '.k','markersize', 30)
hold on
plot(EnergyDensity200, oCAP1_Ampl200, '.r','markersize', 30)
xlabel('Radiant Energy (mJ/cm^2)')
ylabel('uV')
set(gca, 'fontsize',15)
set(gcf, 'Position',[600 500 700 400])
legend('oCAP2','oCAP1','location', 'northwest')
title(strcat(exptdate,' - Mouse 1 - 200 um Fiber'),'interpreter','none')
xlim([0 20])
set(gcf,'Position',[400 400 300 300])

%cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
%exportgraphics(gcf, 'Figure1_PartE_oCAP_PowerCurve_inset.pdf', 'ContentType', 'vector');

figure
plot(EnergyDensity200, oCAP2_Ampl200, 'marker','s','linestyle','none','color',[0.16 0.41 1],'markersize', 12,'linewidth',2)
hold on
plot(EnergyDensity200, oCAP1_Ampl200, 'marker','s','linestyle','none','color',[0 0.7 0],'markersize', 12,'linewidth',2)
plot(EnergyDensity, oCAP2_Ampl, '.','color',[0.16 0.41 1],'markersize', 30)
plot(EnergyDensity, oCAP1_Ampl, '.','color',[0 0.7 0],'markersize', 30)
xlabel('Radiant Energy (mJ/cm^2)')
ylabel('uV')
set(gca, 'fontsize',15)

legend('oCAP2 (200um)', 'oCAP1 (200um)','oCAP2 (400um)','oCAP1(400um)','location', 'southeast')
title(strcat(exptdate,' - Mouse 1 - 400 um & 200 um Fiber'),'interpreter','none')
set(gcf,'Position',[700 800 500 400])
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Supplement_200_vs_400um_fiber.pdf', 'ContentType', 'vector');

% Inset
figure
plot(EnergyDensity200, oCAP2_Ampl200, 'marker','s','linestyle','none','color',[0.16 0.41 1],'markersize', 12,'linewidth',2)
hold on
plot(EnergyDensity200, oCAP1_Ampl200, 'marker','s','linestyle','none','color',[0 0.7 0],'markersize', 12,'linewidth',2)
plot(EnergyDensity, oCAP2_Ampl, '.','color',[0.16 0.41 1],'markersize', 30)
plot(EnergyDensity, oCAP1_Ampl, '.','color',[0 0.7 0],'markersize', 30)
xlabel('Radiant Energy (mJ/cm^2)')
ylabel('uV')
set(gca, 'fontsize',15)
xlim([0 20])

legend('oCAP2 (200um)', 'oCAP1 (200um)','oCAP2 (400um)','oCAP1(400um)','location', 'northwest')
title(strcat(exptdate,' - Mouse 1 - 400 um & 200 um Fiber'),'interpreter','none')
set(gcf,'Position',[700 0 500 400])
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Supplement_200_vs_400um_fiber_Inset.pdf', 'ContentType', 'vector');

%% - Part D - ABR and CAP to laser at 2 powers

close all
exptdate = '2023_3_19';
folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
path = strcat(folder,exptdate,'/');
cd(path)

Fs      = 24410;            % Sampling rate in Hz
tstep   = 1/Fs;
t       = [0:268]; % # of time points
t       = (t.*tstep) + (200e-6);  % Time vector, correct for 200 us offset

% Plot the optical ABR and CAP
clear Pulse_Length
Runs            = [15  20  ];
PowerPercent    = [50 17  ];
Pulse_Length    = [100 100 ];

Target_Distance = 300;
Fiber_Diameter  = 400;

% Make 3-D array containing all the CSV data
for i = 1:length(Runs)
    s(i,:,:) = readmatrix(strcat('Run', num2str(Runs(i)), '.csv'));
end

subplot(1,2,1)
offset = 0;
clear EnergyDensity Watts
for j = 1:length(Runs)
    % Record the power and radiant energy delivered -> only need this once,
    % plot it w/the CAP
    %[Watts(j), ~, ~, EnergyDensity(j)] = CapellaEnergy(Fiber_Diameter, PowerPercent(j), Target_Distance, Pulse_Length(j));
    data = squeeze(s(j,2,49:end));

    plot(t.*1000, 0.001.*(data-offset),'linewidth',2,'color',[0 0 0])
    patch([1 1 1.1 1.1], [-10 5 5 -10],'k','facealpha',0.1)
    hold on
 
    offset = offset + 5000;
    grid on
    xlim([0 11])
    ylim([-10 5])
    set(gca, 'fontsize',13)
end
title('Optical ABR')

subplot(1,2,2)
offset = 0;
clear EnergyDensity
for j = 1:length(Runs)
    % Record the power and radiant energy delivered
    [Watts(j), ~, ~, EnergyDensity(j)] = CapellaEnergy(Fiber_Diameter, PowerPercent(j), Target_Distance, Pulse_Length(j));
    data = squeeze(s(j,3,49:end));

    plot(t.*1000, 0.001.*(data-offset),'linewidth',2,'color',[0 0 0])
    patch([1 1 1.1 1.1], [-100 50 50 -100],'k','facealpha',0.1)
    hold on
 
    offset = offset + 50000;
    grid on
    xlim([0 11])
    ylim([-100 50])     
    set(gca, 'fontsize',13)
end
legend(strcat(num2str(EnergyDensity(1)), 'mJ/cm^2'),...
    strcat(num2str(EnergyDensity(2)), 'mJ/cm^2'),...    
    'location','southeast')

set(gcf, 'Position',[0 0 600 400])
title('Optical CAP 400 um')
xlabel('Time (ms)')

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Figure1_PartD_Optical_ABR_CAP.pdf', 'ContentType', 'vector');


