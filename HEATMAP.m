clear all
close all
clc

%% Load the folder for parsing
myDir = '~/Applications/UTORII_DATA/';
% myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/';
myFiles = dir(fullfile(myDir, '*DI*SOI*.mat'));

%% Parse the folder and calculate CoV for each file

% Initialize variables
IMEP = zeros(length(myFiles),1);
CoV_IMEP = zeros(length(myFiles),1);
DI_quantity = zeros(length(myFiles),1);
DI_duration = zeros(length(myFiles),1);
DI_timing = zeros(length(myFiles),1);
CA50 = zeros(length(myFiles),1);
Ammonia = zeros(length(myFiles),1);

for i = 1:length(myFiles)

    % Load files
    baseName = myFiles(i).name;
    load([myDir, baseName]);

    % Control variables
    DI_quantity(i) = mean(LowSpeed.Pilot_1_Injection_Quantity);
    DI_duration(i) = mean(Cylinder_1_Cycle_Data.Injection_1_Duration);
    DI_timing(i) = -mean(Cylinder_1_Cycle_Data.Injection_1_SOI);
    
    % Dependent variables
    IMEP(i) = mean(Cylinder_1_Cycle_Data.IMEPn);
    CoV_IMEP(i) = std(Cylinder_1_Cycle_Data.IMEPn)/mean(Cylinder_1_Cycle_Data.IMEPn)*100;
    CA50(i) = mean(Cylinder_1_Cycle_Data.CA50);
    Ammonia(i) = mean(LowSpeed.ISB_Fuel_Flow_PFI)/(mean(LowSpeed.ISB_Fuel_Flow+LowSpeed.ISB_Fuel_Flow_PFI))*100;
    fprintf('%s\nCoV IMEP:              %.1f %s\n\n', baseName, std(Cylinder_1_Cycle_Data.IMEPn)/mean(Cylinder_1_Cycle_Data.IMEPn)*100, '%');

end
%% Interpolate values

% Refine mesh
x_fine = min(DI_quantity):0.001:max(DI_quantity);
y_fine = min(DI_timing):0.1:max(DI_timing);
[DI_Q, DI_SOI] = meshgrid(x_fine,y_fine);

% Set up fittype and options
ft = 'cubicinterp';
opts = fitoptions( 'Method', 'CubicSplineInterpolant' );
opts.ExtrapolationMethod = 'none';
opts.Normalize = 'on';

% Fit model to data
CoV_m = fit( [DI_quantity, DI_timing], CoV_IMEP, ft, opts );
CoV_IMEP_spline = CoV_m(DI_Q, DI_SOI);
CA50_m = fit( [DI_quantity, DI_timing], CA50, ft, opts );
CA50_spline = CA50_m(DI_Q, DI_SOI);
Ammonia_m = fit( [DI_quantity, DI_timing], Ammonia, ft, opts );
Ammonia_spline = Ammonia_m(DI_Q, DI_SOI);

%% Plot results

figure(1); clf; hold on
contourf(DI_Q, DI_SOI, Ammonia_spline,85:95,'linestyle','none'); clim([85 95]); colorbar;
plot(DI_quantity, DI_timing,'o','markersize',6,'MarkerFaceColor',0.7*[1 1 1],'MarkerEdgeColor','b');
contour(DI_Q, DI_SOI, CoV_IMEP_spline, [0:3:21 2],':k','linewidth',2,'showtext','on');
contour(DI_Q, DI_SOI, CA50_spline, [8 8],'-k','linewidth',2);
legend('Ammonia content (%)', 'Data points', 'CoV IMEP(%)','CA50 = 8 [CA deg]','location','southeast')
xlabel('Injection Quantity (mg)');
ylabel('Injection timing (deg bTDC)');
hold off