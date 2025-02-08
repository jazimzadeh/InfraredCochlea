function [Watts, Energy, Area, EnergyDensity] = CapellaEnergy(Fiber_Diameter, Percent_Power, Target_Distance, Pulse_Length)
%
% INPUTS
% Fiber_Diameter    - '200' or '400' um optical fiber diameter
% Percent_Power     - % Power chosen on Capella laser
% Target_Distance   - Distance in um between fiber tip and target
% Pulse_Length      - Duration of laser pulse in microseconds
%
% OUTPUTS
% Watts             - Power delivered in Watts (W)
% Area              - Area over which power was delivered (cm^2)
% EnergyDensity     - Energy per area at target distance (mJ/cm^2)

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/MATLAB code/")
load('Capella_Power', 'Power_Commanded', 'Power_400', 'Power_200')

if Fiber_Diameter == 400
    mWatts = interp1(Power_Commanded,Power_400, Percent_Power,'makima');
elseif Fiber_Diameter == 200
    mWatts = interp1(Power_Commanded,Power_200, Percent_Power,'makima');
end

Watts   =  mWatts/1000;

r       = (Target_Distance * tand(7.7)) + Fiber_Diameter/2;     % radius in um
Area    = pi * (r^2);                                           % A in um^2
Area    = Area *10^(-8);                                        % A in cm^2

Energy  = Watts * Pulse_Length * 10^-6;                         % E in Joules
Energy  = Energy * 10^3;                                        % E in mJ
EnergyDensity = Energy/Area;                                    % mJ/cm^2