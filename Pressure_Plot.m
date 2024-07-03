clear all
close all
clc

% Load the file

% filename = "20240305_1200_DI9_SOI42_PFI8_5_020.mat";
% load(filename)
myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/';
filename = '20240305_1200_DI8_4_SOI44_PFI8_5_010.mat';
load([myDir, filename])

% extract the filename for output data
file = erase(filename, ".mat");
ofile = append(file, "_X_res.csv");

% Plot the Cylinder Pressure vs. Crank Angle

Pmax = 11; % MPa
n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);
gamma = 1.3;
EVC = -355;
EVO = 165;
CA_deg = linspace(-359.8, 360, 3600);
Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure, [3600, n_cycles]) * 10 ^ -3;

i_evc = find(CA_deg == -355);
i_evo = find(CA_deg == 160.8000);

figure(1); clf
plot_Pcyl(Pcyl_CA, Pmax, EVC, EVO);

% Calculate the X_res as a function of cycle, k

X_res = zeros([1, 2000]);

for i = 1:2000
    X_res(i) =  (Volume.Volume(i_evc) / Volume.Volume(i_evo)) * (Pcyl_CA(i_evc, i) / Pcyl_CA(i_evo, i)) ^ (1 / gamma);
end
% 2000 cycles, 3600 data points per cycle = 7,200,000
% writematrix(X_res, ofile);

% Function to plot, lifted from ANALYSIS.m
function plot_Pcyl(Pcyl_CA, Pmax, EVC, EVO)
    % Crank angle undisvidual cycle
    CA_deg = linspace(-359.8, 360, 3600);

    % Isolate compression and power stroke
    idx_comb  = find( CA_deg >= -360 & CA_deg <= 360 );
    CAD_comb  = CA_deg(idx_comb);
    
    % Combustion part
    Pcyl_comb = Pcyl_CA(idx_comb,:);
    
    % Plot
    colormap(summer);

    plot(CAD_comb,Pcyl_comb);
    xline(EVC, 'r', {'EVC'});
    xline(EVO, 'r', {'EVO'});
    xlabel('Crank Angle (deg)');
    ylabel('Cylinder Pressure (MPa)'); set(gca,'xtick',-400:60:400);
    xlim([-400 400]); ylim([0 Pmax]);
    title("Cylinder Pressure");
end