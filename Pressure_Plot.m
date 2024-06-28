clear all
close all
clc

% Load the file

load('20240305_1200_DI7_8_SOI43_PFI8_5_019.mat')

% Plot the Cylinder Pressure vs. Crank Angle

Pmax = 11; % MPa
n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);
gamma = 1.3;
EVC = -355;
EVO = 165;

figure(1); clf
plot_Pcyl(Cylinder_1_Synch_Data.Cylinder_Pressure, Pmax, n_cycles, EVC, EVO);

% Calculate the X_res as a function of cycle, k

% Ask Dr. Maldonado/Kaul to explain timing/indexing
X_res = (Volume.Volume() / Volume.Volume()) * (Cylinder_1_Synch_Data.Cylinder_Pressure() / Cylinder_1_Synch_Data.Cylinder_Pressure()) ^ (1 / gamma);

% Function to plot, lifted from ANALYSIS.m
function plot_Pcyl(Pcyl, Pmax, n_cycles, EVC, EVO)
    % Crank angle undividual cycle
    CA_deg = linspace(-359.8, 360, 3600);
    Pcyl_CA = reshape(Pcyl, [3600, n_cycles]) * 1e-3;

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