%% Tabula rasa
clear all
clc

%% Load data
myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/20240305/';
myFiles = dir(fullfile(myDir, '*DI*SOI*.mat'));

% for i = 1:length(myFiles)

i = 15;

% Load files
fileName = myFiles(i).name; 
load([myDir, fileName]);

% Get file name without extension
baseName = erase(fileName, ".mat");
baseName = baseName(15:end-10);

% Estimate fuel quantity

% Number of cycles in file
n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);

% Lower heating values (J)
Q_LHV_diesel  = 44.1 * 1e6;
Q_LHV_ammonia = 18.6 * 1e6;

% % Diesel fuel mass per cycle (Kg)
% DI_quantity = 2 * LowSpeed.ISB_Fuel_Flow ./ LowSpeed.Speed * 1e-3;
% % Get the size of the original signal
% original_size = length(DI_quantity);
% % Create a new index for the resized signal with n_cycles
% new_index = linspace(1, original_size, n_cycles);
% % Use interp1 to interpolate the signal to the new index
% DI_quantity_c = interp1(1:original_size, DI_quantity, new_index)';
% 
% % Introduce diesel injector variability
DI_duration = Cylinder_1_Cycle_Data.Injection_1_Duration';
% DI_quantity_c2c = (DI_duration / median(DI_duration)) .^ (1/2) .* DI_quantity_c;

% Diesel fuel mass per cycle (Kg)
DI_quantity = LowSpeed.Pilot_1_Injection_Quantity * 1e-6;
% Get the size of the original signal
original_size = length(DI_quantity);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
DI_quantity_c2c = interp1(1:original_size, DI_quantity, new_index)';

% Ammonia fuel mass per cycle (Kg)
PI_quantity = 2 * LowSpeed.ISB_Fuel_Flow_PFI ./ LowSpeed.Speed * 1e-3;
% Get the size of the original signal
original_size = length(PI_quantity);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
PI_quantity_c2c = interp1(1:original_size, PI_quantity, new_index)';

% Estimate residual gas fraction

% Model parameters
EVC = -355;
EVO = 165;
gamma = 1.3;

% Crank angle axis
CA_deg = linspace(-359.8, 360, 3600);

% Cylinder pressure (Pa)
if length(Cylinder_1_Synch_Data.Cylinder_Pressure)/3600 == n_cycles
    Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure, [3600, n_cycles]) * 10 ^ -3;
else
    Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure(1:end-3600), [3600, n_cycles]) * 10 ^ -3;
end

% Valve timings
i_evc = find(CA_deg == -355);
i_evo = find(CA_deg == 160.8000);

% Calculate the X_res (-)
X_res = zeros(n_cycles, 1);
for i = 1:n_cycles
    X_res(i) =  (Volume.Volume(i_evc) / Volume.Volume(i_evo)) * (Pcyl_CA(i_evc, i) / Pcyl_CA(i_evo, i)) ^ (1 / gamma);
end

% Estimate combustion efficiency

% Potential heat release
Q_potential = DI_quantity_c2c * Q_LHV_diesel + PI_quantity_c2c * Q_LHV_ammonia;

% Gross heat release (J)
Q_gross = Cylinder_1_Cycle_Data.Gross_Heat_Release';

% Estimated combustion efficiency
eta_c = (1 - X_res) ./ (Q_potential ./ Q_gross - X_res);

% Measured combustion efficiency
eta_c_mea = LowSpeed.Combustion_Efficiency / 100;
% Get the size of the original signal
original_size = length(eta_c_mea);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
eta_c_mea_c2c = interp1(1:original_size, eta_c_mea, new_index);

figure; hold on; box on;
plot(eta_c);
plot(eta_c_mea_c2c, 'LineWidth', 2); legend('Estimated', 'Measured');
ylabel('Combustion efficiency'); xlabel('Cycle')
% print(['Model_Plots/', baseName, 'eta_c'], '-dpng', '-r300');

%% Estimate in-cylinder mass

% Air mass per cycle (Kg)
Air_quantity = 2 * LowSpeed.Mass_Air_Flow ./ LowSpeed.Speed * 1e-3;
% Get the size of the original signal
original_size = length(Air_quantity);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
Air_quantity_c2c = interp1(1:original_size, Air_quantity, new_index)';

% Initialize variables
M_fuel = zeros(n_cycles,1);
M_air = zeros(n_cycles,1);

% LHV for mix
Q_LHV_mix = (DI_quantity_c2c(1) * Q_LHV_diesel + PI_quantity_c2c(1) * Q_LHV_ammonia) ./ (DI_quantity_c2c(1) + PI_quantity_c2c(1));

% Initial condition
M_fuel(1) = ((DI_quantity_c2c(1) + PI_quantity_c2c(1)) - X_res(1) * Q_gross(1) / Q_LHV_mix) / (1 - X_res(1));
M_air(1) = (Air_quantity_c2c(1) + X_res(1) * Q_gross(1) / Q_LHV_mix) / (1 - X_res(1));

% Propagate system forward
for i = 1:n_cycles-1
    state = [M_fuel(i); M_air(i)];
    Matrix_A = X_res(i) * [1 - eta_c(i), 0; eta_c(i), 1];
    input = [DI_quantity_c2c(i) + PI_quantity_c2c(i); Air_quantity_c2c(i)];
    next_state = Matrix_A * state + input;

    M_fuel(i+1) = next_state(1);
    M_air(i+1) = next_state(2);
end

figure
subplot(2,1,1); plot(M_fuel*1e6); ylabel('In-cylinder fuel estimate (mg)')
subplot(2,1,2); plot(M_air*1e6); ylabel('In-cylinder air estimate (mg)')
xlabel('Cycles');
% print(['Model_Plots/', baseName, 'Mass'], '-dpng', '-r300');

%% Residual gas fraction as function of Q_gross

% Residual gas fraction in percentage (%)
X_res_per = X_res * 100;

% Estimate the mean and covariance matrix
X_res_mu    = mean([X_res_per, Q_gross]);
X_res_Sigma = cov([X_res_per, Q_gross]);

% Conditional Gaussian
X_res_per_model = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross);

figure
subplot(2,1,1); hold on; box on
scatter(Q_gross, X_res_per); hold on
scatter(Q_gross, X_res_per_model); legend('Experimental', 'Conditional Gaussian')
xlabel('Q_{Gross} (J)'); ylabel('X_{res} (%)')

subplot(2,1,2); hold on; box on
histogram(X_res_per, "Normalization", "pdf");
histogram(X_res_per_model, "Normalization", "pdf");
legend('Experimental', 'Conditional Gaussian')
xlabel('X_{res} (%)'); ylabel('PDF')

% print(['Model_Plots/', baseName, 'X_res_model'], '-dpng', '-r300');

%% Combustion Efficiency as function of AFR

% Combustion efficiency in percentage (%)
eta_c_per = eta_c * 100;

% State
Mass_mg = [M_fuel, M_air] * 1e6; 

% Estimate the mean and covariance matrix
eta_c_mu    = mean([eta_c_per, Mass_mg]);
eta_c_Sigma = cov([eta_c_per,Mass_mg]);

% Conditional Gaussian
eta_c_per_model = conditional_Gauss(eta_c_mu, eta_c_Sigma, Mass_mg);

figure
subplot(2,1,1); hold on; box on
scatter(M_air./M_fuel, eta_c_per); ylim([0, 100]);
scatter(M_air./M_fuel, eta_c_per_model); legend('Experimental', 'Conditional Gaussian')
xlabel('Air-to-Fuel Ratio (-)'); ylabel('Combustion Efficiency (%)')

subplot(2,1,2); hold on; box on
histogram(eta_c_per, "Normalization", "pdf")
histogram(eta_c_per_model, "Normalization", "pdf");
legend('Experimental', 'Conditional Gaussian')
xlabel('Combustion Efficiency (%)'); ylabel('PDF')

% print(['Model_Plots/', baseName, 'eta_c_model'], '-dpng', '-r300');

%% Simulate dynamic system

% Initialize variables
M_fuel_sim = zeros(n_cycles,1);
M_air_sim = zeros(n_cycles,1);
eta_c_sim = zeros(n_cycles,1);
X_res_sim = zeros(n_cycles,1);
Q_gross_sim = zeros(n_cycles,1);

% Initial condition
M_fuel_sim(1) = M_fuel(1);
M_air_sim(1) = M_air(1);

% Inputs
Diesel_fuel = DI_quantity_c2c;
Ammonia_fuel = PI_quantity_c2c;
Fresh_air = Air_quantity_c2c;

for i = 1:n_cycles
    
    % State
    state = [M_fuel_sim(i); M_air_sim(i)];
    
    % Combustion efficiency
    eta_c_sim(i) = conditional_Gauss(eta_c_mu, eta_c_Sigma, state'*1e6) / 100;
    
    % Effective LHV
    Q_LHV_eff = (Diesel_fuel(i) * Q_LHV_diesel + Ammonia_fuel(i) * Q_LHV_ammonia) ./ (Diesel_fuel(i) + Ammonia_fuel(i));

    % Gross heat release
    Q_gross_sim(i) = eta_c_sim(i) * M_fuel_sim(i) * Q_LHV_eff;

    % Residual gas fraction
    X_res_sim(i) = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim(i)) / 100;

    % Residual mass matrix
    Matrix_res = X_res_sim(i) * [1 - eta_c_sim(i), 0; eta_c_sim(i), 1];

    % Fresh fuel and air
    input = [Diesel_fuel(i) + Ammonia_fuel(i); Fresh_air(i)];

    if i < n_cycles
        % Calculate next cycle
        next_state = Matrix_A * state + input;
        M_fuel_sim(i+1) = next_state(1);
        M_air_sim(i+1) = next_state(2);
    end
end

% Fix initial condition
M_fuel_sim(1) = M_fuel_sim(2);
M_air_sim(1) = M_air_sim(2);

figure
subplot(4,2,1); hold on; box on;
plot(M_fuel*1e6); plot(M_fuel_sim*1e6); ylabel('M_{fuel} (mg)');
xLimits = ylim; title('Timeseries')
subplot(4,2,2); hold on; box on; xlim(xLimits)
histogram(M_fuel*1e6, "Normalization", "pdf")
histogram(M_fuel_sim*1e6, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse'); title('PDF')

subplot(4,2,3); hold on; box on;
plot(M_air*1e6); plot(M_air_sim*1e6); ylabel('M_{air} (mg)'); xLimits = ylim;
subplot(4,2,4); hold on; box on; xlim(xLimits)
histogram(M_air*1e6, "Normalization", "pdf")
histogram(M_air_sim*1e6, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

subplot(4,2,5); hold on; box on;
plot(Q_gross); plot(Q_gross_sim); ylabel('Q_{gross} (J)'); xLimits = ylim;
subplot(4,2,6); hold on; box on; xlim(xLimits)
histogram(Q_gross, "Normalization", "pdf")
histogram(Q_gross_sim, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

subplot(4,2,7); hold on; box on;
plot(X_res*100); plot(X_res_sim*100); ylabel('X_{res} (%)'); xLimits = ylim;
subplot(4,2,8); hold on; box on; xlim(xLimits)
histogram(X_res*100, "Normalization", "pdf")
histogram(X_res_sim*100, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

% print(['Model_Plots/', baseName, 'Simulator'], '-dpng', '-r300');

%% Combustion phasing

CA50 = Cylinder_1_Cycle_Data.CA50';
DI_timing = Cylinder_1_Cycle_Data.Injection_1_SOI';

%% Desired operating conditions

DI_quantity_des = mean(LowSpeed.Pilot_1_Injection_Quantity);
DI_timing_des = -mean(Cylinder_1_Cycle_Data.Injection_1_SOI);

%% Safe file

% model_data = table(M_fuel_sim, M_air_sim, eta_c, Q_gross, X_res, Diesel_fuel, ...
%     Ammonia_fuel, Fresh_air, CA50, -DI_timing, DI_duration);
% writetable(model_data, ['Model_Data/', baseName, 'Data.csv']);
% save(['Model_Data/', baseName, 'Parameters.mat'], 'Q_LHV_diesel', ...
%     'Q_LHV_ammonia', 'eta_c_mu', 'eta_c_Sigma', 'X_res_mu', 'X_res_Sigma', ...
%     'DI_quantity_des', 'DI_timing_des');

% end
%% Auxiliary functions

function x1_sim = conditional_Gauss(mu, Sigma, x2)

% Conditional Gaussian
% x1 | x2 ~ N(mu_cond, sigma_cond)

% Conditional standard deviation
sigma_cond = sqrt(Sigma(1, 1) - Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * Sigma(2:end, 1));

% Gaussian Conditional mean
mu_cond = mu(1) + Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * (x2 - mu(2:end))';

% Simulate residual gas fraction in percentage
x1_sim = normrnd(mu_cond, sigma_cond);
end