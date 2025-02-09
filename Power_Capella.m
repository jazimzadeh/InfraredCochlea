%% - Import the measured power data
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets/")
Power = readmatrix(strcat('2023_3_8_Capella_Power_Measurements.csv'));

%% - Correct for duty cycle 
% - Data was measured with laser on 5% of the time
% - So divide by 0.05 to get true peak power value
Power_400       = Power(:,2)./0.05;
Power_200       = Power(:,3)./0.05;
Power_Aperture  = Power(:,4)./0.05;
Power_Commanded = Power(:,1);

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/MATLAB code/")
save('Capella_Power') % Save this data to a MAT file in my MATLAB code folder

%% - Plot Power vs. Power commanded (%)
clf
plot(Power(:,1), Power_400, '.','markersize',30)
hold on
plot(Power(:,1), Power_200, '.','markersize',30)
plot(Power(:,1), Power_Aperture, '.','markersize',30)
plot(Power(:,1), Power_200.*2.918, 'o','markersize',11,'linewidth',2)
plot(Power(:,1), Power_Aperture.*0.92, 'or','markersize',11,'linewidth',2)
set(gca,'fontsize',16)
set(gcf,'Position',[0 100 900 800])
ylabel('Power (mW)')
xlabel('Power Commanded (%)')
title('Power (mW) vs Power Commanded (%)')

legend('400 um Fiber', ...
    '200 um Fiber', ...
    'Capella Aperture', ...
    '200 um Fiber x 2.918', ...
    'Capella Aperture x 0.92', ...
    'location','northwest', ...
    'fontsize',16)

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, '_Beam_Power_vs_CommandedPower.pdf', 'ContentType', 'vector');

%% - Interpolation
% Use Akima method of interpolation: 
% [1] Akima, Hiroshi. "A new method of interpolation and smooth curve fitting based on local procedures." Journal of the ACM (JACM) , 17.4, 1970, pp. 589-602.
% [2] Akima, Hiroshi. "A method of bivariate interpolation and smooth surface fitting based on local procedures." Communications of the ACM , 17.1, 1974, pp. 18-20.

Interp_Powers = (1:0.5:100); % Vector of commanded power % for which to generate interpolated data
Interp_400 = interp1(Power_Commanded,Power_400, Interp_Powers ,'makima'); % Akima interpolation 
Interp_400l = interp1(Power_Commanded,Power_400, Interp_Powers ,'linear'); % Linear interpolation 

clf
plot(Power_Commanded, Power_400, '.', 'markersize', 40) % Plot original data
hold on
plot(Interp_Powers, Interp_400, '.', 'markersize', 14) % Plot Akima interpolated data
plot(Interp_Powers, Interp_400l, '.', 'markersize', 10) % Plot Linear interpolated data
title('Linear vs Akima Interpolation')
ylabel('Power (mW)')
xlabel('Power Commanded (%)')
set(gca,'fontsize',16)

legend('400 um Fiber Data', ...
    'Akima Interpolated Data (by 0.5%)',...
    'Linear Interpolated Data (by 0.5%)',...
    'location','northwest', ...
    'fontsize',16)

%% - Usage Example: Output energy, given % Power
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/MATLAB code/")
load('Capella_Power')

% Usage example, when power at a particular % is desired:
Laser_Command = 50.15;
Interp_Rand= interp1(Power_Commanded,Power_400, Laser_Command,'makima');

%% - Import the measured beam width data - 400 um & 200 um fibers
cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets/")
close all
Beam400 = readmatrix(strcat('2023_3_8_Capella_400um_Beam.csv'));
Beam200 = readmatrix(strcat('2023_3_8_Capella_200um_Beam.csv'));

Razor_pos_400 = Beam400(:,1) - min(Beam400(:,1));
Razor_pos_200 = Beam200(:,1) - min(Beam200(:,1));

% Correct for duty cycle and interpolate - 400 um fiber
Interp_Pos400_0mm   = (900:2:1500); % Vector of commanded positions for which to generate interpolated data
Interp_Power400_0mm = interp1(Razor_pos_400, Beam400(:,2)./0.05, Interp_Pos400_0mm ,'makima'); % Akima interpolation 
Interp_Pos400_2mm   = (0:5:2400); % Vector of commanded positions for which to generate interpolated data
Interp_Power400_2mm = interp1(Razor_pos_400, Beam400(:,3)./0.05, Interp_Pos400_2mm ,'makima'); % Akima interpolation 

% Correct for duty cycle and interpolate - 200 um fiber
Interp_Pos200_0mm   = (340:1:600); % Vector of commanded positions for which to generate interpolated data
Interp_Power200_0mm = interp1(Razor_pos_200, Beam200(:,2)./0.05, Interp_Pos200_0mm ,'makima'); % Akima interpolation 
Interp_Pos200_2mm   = (0:4:1300); % Vector of commanded positions for which to generate interpolated data
Interp_Power200_2mm = interp1(Razor_pos_200, Beam200(:,3)./0.05, Interp_Pos200_2mm ,'makima'); % Akima interpolation 

% Plot original data, corrected for duty cycle - 400 um Fiber
clf
plot(Razor_pos_400, Beam400(:,2)./0.05,'.', 'markersize', 30)
hold on
plot(Razor_pos_400, Beam400(:,3)./0.05,'.', 'markersize', 30)
% Overlay interpolated data
plot(Interp_Pos400_0mm, Interp_Power400_0mm, '.')
plot(Interp_Pos400_2mm, Interp_Power400_2mm, '.')
title('400um Fiber - Beam width and spread')
ylabel('Power (mW)')
xlabel('Razorblade Position (um)')
set(gca,'fontsize',16)
set(gcf, 'Position', [0 300 700 600])

legend('400 um Fiber, 0 mm to razor', ...
    '400 um Fiber, 2 mm to razor',...
    'Interpolation',...
    'Interpolation',...
    'location','northeast', ...
    'fontsize',16)

% Plot original data, corrected for duty cycle - 200 um Fiber
figure
plot(Razor_pos_200, Beam200(:,2)./0.05,'.', 'markersize', 30)
hold on
plot(Razor_pos_200, Beam200(:,3)./0.05,'.', 'markersize', 30)
% Overlay interpolated data
plot(Interp_Pos200_0mm, Interp_Power200_0mm, '.')
plot(Interp_Pos200_2mm, Interp_Power200_2mm, '.')
title('200um Fiber - Beam width and spread')
ylabel('Power (mW)')
xlabel('Razorblade Position (um)')
set(gca,'fontsize',16)
ylim([0 1600])
set(gcf, 'Position', [700 300 700 600])

legend('200 um Fiber, 0 mm to razor', ...
    '200 um Fiber, 2 mm to razor',...
    'Interpolation',...
    'Interpolation',...
    'location','northeast', ...
    'fontsize',16)

% - Flip the data so power starts at 0
Pos400_2mm = -(Interp_Pos400_2mm - max(Interp_Pos400_2mm));
Pos400_0mm = -(Interp_Pos400_0mm - max(Interp_Pos400_2mm));

Pos200_2mm = -(Interp_Pos200_2mm - max(Interp_Pos200_2mm));
Pos200_0mm = -(Interp_Pos200_0mm - max(Interp_Pos200_2mm));

figure
plot(Pos400_2mm, Interp_Power400_2mm, '.')
hold on
plot(Pos400_0mm, Interp_Power400_0mm, '.')
set(gcf,'Position',[800 0 500 300])

%% - Smooth and find derivative to get beam profile - 400 um Fiber
n_smooth = 50;
Power400_2mm_smooth = smooth(Interp_Power400_2mm, n_smooth);
Power400_0mm_smooth = smooth(Interp_Power400_0mm, n_smooth);

dydx = diff(Power400_2mm_smooth)./diff(Pos400_2mm);
halfmax = max(dydx)/2;
dydx_0 = diff(Power400_0mm_smooth)./diff(Pos400_0mm);
halfmax_0 = max(dydx_0)/2;

close all
subplot(2,1,1)
plot(Pos400_2mm, Interp_Power400_2mm, '.k')
hold on
plot(Pos400_2mm, Power400_2mm_smooth, '.')
plot(Pos400_0mm, Interp_Power400_0mm, '.k')
plot(Pos400_0mm, Power400_0mm_smooth, '.')

legend('2 mm, Original',...
    '2 mm, Smoothed',...
    '0 mm, Original',...
    '0 mm, Smoothed',...
    'location','northwest')
title('400 um Fiber')
set(gca, 'fontsize', 16)
ylabel('Power (mW)')

subplot(2,1,2)
plot(Pos400_2mm(1:end-1), dydx(:,1), '.k') % Plot derivative
hold on
plot(Pos400_0mm(1:end-1), dydx_0(:,1), '.k') % Plot derivative
ylim([0 8])
set(gca, 'fontsize', 16)

[idxMax, ~]     = IntersectionIndex(dydx(1:end/2,1), halfmax(1,1));
[idxMax_0, ~]   = IntersectionIndex(dydx_0(1:end/2,1), halfmax_0(1,1));
plot(Pos400_2mm(idxMax), dydx(idxMax,1), 'or','markersize',15,'linewidth',2)
plot(Pos400_0mm(idxMax_0), dydx_0(idxMax_0,1), 'or','markersize',15,'linewidth',2)

[idxMin, ~]     = IntersectionIndex(dydx(end/2:end,1), halfmax(1,1));
[idxMin_0, ~]   = IntersectionIndex(dydx_0(end/2:end,1), halfmax_0(1,1));
idxMin          = idxMin + length(dydx)/2;
idxMin_0        = idxMin_0 + length(dydx_0)/2;
plot(Pos400_2mm(idxMin), dydx(idxMin,1), 'or','markersize',15,'linewidth',2)
plot(Pos400_0mm(idxMin_0), dydx_0(idxMin_0,1), 'or','markersize',15,'linewidth',2)
xlabel('Position (um)')
ylabel('d(Power) / d(Position)')

FWHM_400_2mm = Pos400_2mm(idxMax) - Pos400_2mm(idxMin);
FWHM_400_0mm = Pos400_0mm(idxMax_0) - Pos400_0mm(idxMin_0);
text(100,6,['FWHM (2 mm) =  ', num2str(FWHM_400_2mm), ' um'],'fontsize',16)
text(100,7,['FWHM (0 mm) =  ', num2str(FWHM_400_0mm), ' um'],'fontsize',16)

divergence = atand(((FWHM_400_2mm/2)-(FWHM_400_0mm/2)) / 2000);


title(strcat('Divergence (degrees) =  ', num2str(divergence)))
set(gcf, 'Position', [0 300 750 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, '_BeamWidth_400um.pdf', 'ContentType', 'vector');

%% - Smooth and find derivative to get beam profile - 200 um Fiber
n_smooth = 90;
Power200_2mm_smooth = smooth(Interp_Power200_2mm, n_smooth);
Power200_0mm_smooth = smooth(Interp_Power200_0mm, n_smooth);

dydx = diff(Power200_2mm_smooth)./diff(Pos200_2mm);
halfmax = max(dydx)/2;
dydx_0 = diff(Power200_0mm_smooth)./diff(Pos200_0mm);
halfmax_0 = max(dydx_0)/2;

figure
subplot(2,1,1)
plot(Pos200_2mm, Interp_Power200_2mm, '.k')
hold on
plot(Pos200_2mm, Power200_2mm_smooth, '.')
plot(Pos200_0mm, Interp_Power200_0mm, '.k')
plot(Pos200_0mm, Power200_0mm_smooth, '.')

legend('2 mm, Original',...
    '2 mm, Smoothed',...
    '0 mm, Original',...
    '0 mm, Smoothed',...
    'location','northwest')
title('200 um Fiber')
set(gca, 'fontsize', 16)
ylabel('Power (mW)')

subplot(2,1,2)
plot(Pos200_2mm(1:end-1), dydx(:,1), '.k') % Plot derivative
hold on
plot(Pos200_0mm(1:end-1), dydx_0(:,1), '.k') % Plot derivative
ylim([0 12])
set(gca, 'fontsize', 16)

[idxMax, ~]     = IntersectionIndex(dydx(1:end/2,1), halfmax(1,1));
[idxMax_0, ~]   = IntersectionIndex(dydx_0(1:end/2,1), halfmax_0(1,1));
plot(Pos200_2mm(idxMax), dydx(idxMax,1), 'or','markersize',15,'linewidth',2)
plot(Pos200_0mm(idxMax_0), dydx_0(idxMax_0,1), 'or','markersize',15,'linewidth',2)

[idxMin, ~]     = IntersectionIndex(dydx(end/2:end,1), halfmax(1,1));
[idxMin_0, ~]   = IntersectionIndex(dydx_0(end/2:end,1), halfmax_0(1,1));
idxMin          = idxMin + round(length(dydx)/2);
idxMin_0        = idxMin_0 + round(length(dydx_0)/2);
plot(Pos200_2mm(idxMin), dydx(idxMin,1), 'or','markersize',15,'linewidth',2)
plot(Pos200_0mm(idxMin_0), dydx_0(idxMin_0,1), 'or','markersize',15,'linewidth',2)
xlabel('Position (um)')
ylabel('d(Power) / d(Position)')

FWHM_200_2mm = Pos200_2mm(idxMax) - Pos200_2mm(idxMin);
FWHM_200_0mm = Pos200_0mm(idxMax_0) - Pos200_0mm(idxMin_0);
text(100,6,['FWHM (2 mm) =  ', num2str(FWHM_200_2mm), ' um'],'fontsize',16)
text(100,7,['FWHM (0 mm) =  ', num2str(FWHM_200_0mm), ' um'],'fontsize',16)

divergence = atand(((FWHM_200_2mm/2)-(FWHM_200_0mm/2)) / 2000);

title(strcat('Divergence (degrees) =  ', num2str(divergence)))
set(gcf, 'Position', [750 300 750 800])

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, '_BeamWidth_200um.pdf', 'ContentType', 'vector');

%% - Compare energy density of both fibers as function of distance

Distance = [0:20:500]; % In um
PowerPercent    = 100; 
PulseLength     = 100; % us

for j = 1:length(Distance)
    [~, ~, ~, EnergyDensity_400(j)] = CapellaEnergy(400, PowerPercent, Distance(j), PulseLength);
    [~, ~, ~, EnergyDensity_200(j)] = CapellaEnergy(200, PowerPercent, Distance(j), PulseLength);
end

%[Watts, Energy, Area, EnergyDensity] = CapellaEnergy(400, 100, 0, 100)
%[Watts, Energy, Area, EnergyDensity] = CapellaEnergy(200, 100, 0, 100)
close all
plot(Distance, EnergyDensity_400, '.k', 'markersize',30)
hold on
plot(Distance, EnergyDensity_200, '.r', 'markersize', 30)
plot(300, 264,'ok','markersize',30,'linewidth',3)
text(280,290,'Tan 2018','fontsize',14)

set(gca,'fontsize',16)
set(gcf,'Position',[0 300 700 600])
xlabel('Distance from fiber tip (um)')
ylabel('Energy Density (mJ/cm^2)')
title('Energy Density Comparison: 200 um and 400 um Fibers')
legend('400 um Fiber', '200 um Fiber')

%% - Comparison to Tan 2018
% They used a 200 um fiber at max power, and got 164 uJ of energy out
% Assuming 300 um target distance, their irradiation area is 6.2e-4 cm^2:
clc
[~, ~, Area, ~] = CapellaEnergy(200, 100, 300, 100)
% This would yield an energy density of 0.164 mJ / (6.2e-4 cm^2) = 
0.164 / (6.2e-4)
% = 264 mJ/cm^2

% [Watts, Energy, Area, EnergyDensity] = CapellaEnergy(Fiber_Diameter, Percent_Power, Target_Distance, Pulse_Length)
%[Watts, Energy, Area, EnergyDensity] = CapellaEnergy(200, 100, 300, 100)
