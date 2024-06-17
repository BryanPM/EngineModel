clear all
close all
clc

load('20240305_1200_DI8_8_SOI43_PFI8_5_014.mat')

%% Print conditions

fprintf('\nLow-speed variables:\n')
fprintf('----------------------\n')
fprintf('Engine speed:          %.0f %s\n', mean(LowSpeed.Speed), LowSpeed.Properties.Speed.Units);
fprintf('Brake Power:           %.1f %s\n', mean(LowSpeed.Energy_Flux_Brake_Power), LowSpeed.Properties.Energy_Flux_Brake_Power.Units);
fprintf('Torque:                %.1f %s\n', mean(LowSpeed.Torque), LowSpeed.Properties.Torque.Units);
fprintf('Pre-Cat Lambda:        %.2f %s\n', mean(LowSpeed.Exhaust_Pre_Cat_Lambda_UEGO), LowSpeed.Properties.Exhaust_Pre_Cat_Lambda_UEGO.Units);
fprintf('Post-Cat Lambda:       %.2f %s\n', mean(LowSpeed.Exhaust_Post_Cat_Lambda_UEGO), LowSpeed.Properties.Exhaust_Post_Cat_Lambda_UEGO.Units);
fprintf('Diesel fuel flow:      %.2f %s\n', mean(LowSpeed.ISB_Fuel_Flow), LowSpeed.Properties.ISB_Fuel_Flow.Units);
fprintf('Ammonia fuel flow:     %.2f %s\n', mean(LowSpeed.ISB_Fuel_Flow_PFI), LowSpeed.Properties.ISB_Fuel_Flow_PFI.Units);
fprintf('Ammonia content:       %.1f %s\n', mean(LowSpeed.ISB_Fuel_Flow_PFI)/(mean(LowSpeed.ISB_Fuel_Flow+LowSpeed.ISB_Fuel_Flow_PFI))*100, '%');
fprintf('Coolant temperature:   %.1f %s\n', mean(LowSpeed.ISB_T_Coolant_HX_In), LowSpeed.Properties.ISB_T_Coolant_HX_In.Units);
fprintf('Oil temperature:       %.1f %s\n', mean(LowSpeed.ISB_T_Oil_HX_In), LowSpeed.Properties.ISB_T_Oil_HX_In.Units);
fprintf('Intake temperature:    %.1f %s\n', mean(LowSpeed.ISB_T_Int_Surge_Tank), LowSpeed.Properties.ISB_T_Int_Surge_Tank.Units);
fprintf('Manifold Pressure:     %.1f %s\n', mean(LowSpeed.MAP), LowSpeed.Properties.MAP.Units);
fprintf('Air mass flow:         %.1f %s\n', mean(LowSpeed.Mass_Air_Flow), LowSpeed.Properties.Mass_Air_Flow.Units);
fprintf('PFI duration:          %.1f %s\n', mean(LowSpeed.PFI_Duration), LowSpeed.Properties.PFI_Duration.Units);
fprintf('PFI timing:            %.0f %s\n', mean(LowSpeed.PFI_Timing), LowSpeed.Properties.PFI_Timing.Units);
fprintf('DI quantity:           %.1f %s\n', mean(LowSpeed.Pilot_1_Injection_Quantity), LowSpeed.Properties.Pilot_1_Injection_Quantity.Units);
fprintf('DI timing:             %.0f %s\n', mean(LowSpeed.Pilot_1_Injection_SOI), LowSpeed.Properties.Pilot_1_Injection_SOI.Units);

fprintf('\nCycle-based variables:\n')
fprintf('----------------------\n')
fprintf('CA50:                  %.1f %s\n', mean(Cylinder_1_Cycle_Data.CA50), Cylinder_1_Cycle_Data.Properties.CA50.Units);
fprintf('Gross Heat Release:    %.1f %s\n', mean(Cylinder_1_Cycle_Data.Gross_Heat_Release), Cylinder_1_Cycle_Data.Properties.Gross_Heat_Release.Units);
fprintf('Net IMEP:              %.1f %s\n', mean(Cylinder_1_Cycle_Data.IMEPn), Cylinder_1_Cycle_Data.Properties.IMEPn.Units);
fprintf('CoV IMEP:              %.1f %s\n', std(Cylinder_1_Cycle_Data.IMEPn)/mean(Cylinder_1_Cycle_Data.IMEPn)*100, '%');
fprintf('DI duration:           %.1f %s\n', mean(Cylinder_1_Cycle_Data.Injection_1_Duration), Cylinder_1_Cycle_Data.Properties.Injection_1_Duration.Units);

%% Cycle to cycle plots

figure(1); clf
plot_histrogram(Cylinder_1_Cycle_Data.CA50, 'CA50')
figure(2); clf
plot_histrogram(Cylinder_1_Cycle_Data.Gross_Heat_Release, 'Gross Heat Release')
figure(3); clf
plot_histrogram(Cylinder_1_Cycle_Data.IMEPn, 'Net IMEP')

figure(4); clf
plot_return_map(Cylinder_1_Cycle_Data.Gross_Heat_Release, 'Gross Heat Release')
figure(5); clf
plot_return_map(Cylinder_1_Cycle_Data.IMEPn, 'Net IMEP')

%% Cycle to cycle plots

% Number of cycles
n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);
Pmax = 15; % Mpa

figure(6); clf
plot_Pcyl(Cylinder_1_Synch_Data.Cylinder_Pressure, Pmax, n_cycles)

figure(7); clf
plot_CA_signal(Cylinder_1_Synch_Data.Mass_Fraction_Burned, n_cycles, ...
    [-30, 60], 'MFB')

%%
function plot_Pcyl(Pcyl, Pmax, n_cycles)
    % Crank angle undividual cycle
    CA_deg = linspace(-359.8, 360, 3600);
    Pcyl_CA = reshape(Pcyl, [3600, n_cycles]) * 1e-3;

    % Isolate compression and power stroke
    idx_comb  = find( CA_deg >= -180 & CA_deg <= 180 );
    CAD_comb  = CA_deg(idx_comb);
    
    % Combustion part
    Pcyl_comb = Pcyl_CA(idx_comb,:);
    % Discretize cylinder pressure
    d_Pcyl = 0:0.1:Pmax;
    % Order pairs
    x_CAD  = kron( ones(size(Pcyl_comb,2),1), CAD_comb' );
    x_Pcyl = Pcyl_comb(:);
    % Histogram
    N_cnt = histcounts2(x_CAD,x_Pcyl,[ numel(idx_comb), numel(d_Pcyl) ], ...
       'YBinLimits', [ min(d_Pcyl) max(d_Pcyl) ], 'Normalization', 'pdf' );
    % Axis
    [X_CAD, X_Pcyl] = meshgrid( CAD_comb, d_Pcyl);
    % Representative cycle
    Pcyl_m = mean(Pcyl_comb,2);
    [~,idx_Pcyl_rep] = min( sum((Pcyl_comb - Pcyl_m).^2) );
    
    % Plot
    colormap(summer);
    set(gcf,'units','centimeters','position',[0,0,13.4,18]);
    % All cycles
    ah1 = subplot(2,1,1); 
    plot(CAD_comb,Pcyl_comb);
    ylabel('Cylinder Pressure (MPa)'); set(gca,'xtick',-180:60:180);
    xlim([CAD_comb(1) CAD_comb(end)]); ylim([0 Pmax]);
    title("Cylinder Pressure");
    % PDF
    ah2 = subplot(2,1,2); hold on; box on;
    contourf(X_CAD,X_Pcyl,N_cnt',100,'linestyle','none')
    p1=plot(CAD_comb,Pcyl_comb(:,idx_Pcyl_rep),'-k','linewidth',1.5);
    xlim([CAD_comb(1) CAD_comb(end)]); ylim([0 Pmax]);
    hcb=colorbar; hcb.Title.String = "PDF"; set(gca,'xtick',-180:60:180);
    set(gca,'ColorScale','log'); legend(p1,'Representative cycle')
    xlabel('Crank Angle (deg)'); ylabel('Cylinder Pressure (MPa)');
    % Adjust position
    pos1 = get(ah1,'Position'); pos2 = get(ah2,'Position');
    pos1(3) = pos2(3); pos1(2) = pos1(2)-0.06; set(ah1,'Position',pos1)
end

function plot_CA_signal(data, n_cycles, x_lim, name)
    CA_deg = linspace(-359.8, 360, 3600);
    data_CA = reshape(data, [3600, n_cycles]);
    plot(CA_deg, data_CA);
    xlim(x_lim)
    xlabel('Crank Angle (deg)');
    ylabel(name);
end

function plot_histrogram(data, name)
    pd = fitdist(data', 'Normal');
    histogram(data,'LineStyle','none','Normalization','pdf');
    hold on;
    x = linspace(min(data), max(data));
    y = pdf(pd, x);
    plot(x, y, 'r', 'LineWidth', 2);
    legend('Experiment', 'Fitted Gaussian');
    xlabel(name);
    ylabel('Probability Density');
end

function plot_return_map(data, name)
    % data_norm = data ./ median(data);
    data_norm = data;
    scatter(data_norm(1:end-1), data_norm(2:end), 'filled');
    xlabel([name, ' [k]']); ylabel([name, ' [k+1]']);
    axis equal; box on;
end
