%% Tabula rasa
clear all
close all
clc

%% Load data
myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/';
myFiles = dir(fullfile(myDir, '*DI*SOI*.mat'));

% for i = 1:length(myFiles)
i=24;

% Load files
baseName = myFiles(i).name;
load([myDir, baseName]);

%% Estimate fuel quantity

% Number of cycles in file
n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);

% Lower heating values (J)
Q_LHV_diesel  = 44.1 * 1e6;
Q_LHV_ammonia = 18.6 * 1e6;

% Diesel fuel mass per cycle (Kg)
DI_quantity = 2 * LowSpeed.ISB_Fuel_Flow ./ LowSpeed.Speed * 1e-3;
% Get the size of the original signal
original_size = length(DI_quantity);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
DI_quantity_c = interp1(1:original_size, DI_quantity, new_index);

% Introduce diesel injector variability
DI_duration = Cylinder_1_Cycle_Data.Injection_1_Duration;
DI_quantity_c2c = (DI_duration / median(DI_duration)) .^ (1/2) .* DI_quantity_c;

% Ammonia fuel mass per cycle (Kg)
PI_quantity = 2 * LowSpeed.ISB_Fuel_Flow_PFI ./ LowSpeed.Speed * 1e-3;
% Get the size of the original signal
original_size = length(PI_quantity);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
PI_quantity_c2c = interp1(1:original_size, PI_quantity, new_index);

%% Estimate residual gas fraction

% Model parameters
EVC = -355;
EVO = 165;
gamma = 1.3;

% Crank angle axis
CA_deg = linspace(-359.8, 360, 3600);

% Cylinder pressure (Pa)
Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure, [3600, n_cycles]) * 10 ^ -3;

% Valve timings
i_evc = find(CA_deg == -355);
i_evo = find(CA_deg == 160.8000);

% Calculate the X_res (-)
X_res = zeros([1, n_cycles]);
for i = 1:n_cycles
    X_res(i) =  (Volume.Volume(i_evc) / Volume.Volume(i_evo)) * (Pcyl_CA(i_evc, i) / Pcyl_CA(i_evo, i)) ^ (1 / gamma);
end

%% Estimate combustion efficiency

% Estimate combined Q_LHV
Q_LHV_mix = (DI_quantity_c2c * Q_LHV_diesel + PI_quantity_c2c * Q_LHV_ammonia) ./ (DI_quantity_c2c + PI_quantity_c2c);

% Gross heat release (J)
Q_gross = Cylinder_1_Cycle_Data.Gross_Heat_Release;

% Estimated combustion efficiency
eta_c = (1 - X_res) ./ ((DI_quantity_c2c + PI_quantity_c2c) .* Q_LHV_mix ./ Q_gross - X_res);

% Measured combustion efficiency
eta_c_mea = LowSpeed.Combustion_Efficiency / 100;
% Get the size of the original signal
original_size = length(eta_c_mea);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
eta_c_mea_c2c = interp1(1:original_size, eta_c_mea, new_index);

figure; hold on
plot(eta_c)
plot(eta_c_mea_c2c, 'LineWidth', 2)
ylabel('Combustion efficiency')
legend('Estimated', 'Measured')

%% Estimate in-cylinder mass

% Air mass per cycle (Kg)
Air_quantity = 2 * LowSpeed.Mass_Air_Flow ./ LowSpeed.Speed * 1e-3;
% Get the size of the original signal
original_size = length(Air_quantity);
% Create a new index for the resized signal with n_cycles
new_index = linspace(1, original_size, n_cycles);
% Use interp1 to interpolate the signal to the new index
Air_quantity_c2c = interp1(1:original_size, Air_quantity, new_index);

% Initialize variables
M_fuel = zeros(n_cycles,1);
M_air = zeros(n_cycles,1);
Q_gross_model = zeros(n_cycles,1);

% Initial condition
M_fuel(1) = ((DI_quantity_c2c(1) + PI_quantity_c2c(1)) - X_res(1) * Q_gross(1) / Q_LHV_mix(1)) / (1 - X_res(1));
M_air(1) = (Air_quantity_c2c(1) + X_res(1) * Q_gross(1) / Q_LHV_mix(1)) / (1 - X_res(1));
Q_gross_model(1) = eta_c(1) * M_fuel(1) * Q_LHV_mix(1);

% Propagate system forward
for i = 1:n_cycles-1
    state = [M_fuel(i); M_air(i)];
    Matrix_A = X_res(i) * [1 - eta_c(i), 0; eta_c(i), 1];
    input = [DI_quantity_c2c(i) + PI_quantity_c2c(i); Air_quantity_c2c(i)];
    next_state = Matrix_A * state + input;

    M_fuel(i+1) = next_state(1);
    M_air(i+1) = next_state(2);
    Q_gross_model(i+1) = eta_c(i+1) * M_fuel(i+1) * Q_LHV_mix(i+1);
end

figure; plot(M_fuel)
figure; plot(M_air)
figure; hold on; plot(Q_gross_model), plot(Q_gross)

%% Residual gas fraction as function of Q_gross

% Residual gas fraction in percentage (%)
X_res_per = X_res' * 100;

% Estimate the mean and covariance matrix
X_res_mu    = mean([X_res_per, Q_gross']);
X_res_Sigma = cov([X_res_per, Q_gross']);

% Conditional Gaussian
X_res_per_sim = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross');

figure; hold on
scatter(Q_gross, X_res_per); hold on
scatter(Q_gross, X_res_per_sim); legend('experiment', 'simulation')
xlabel('Q_{Gross} (J)'); ylabel('X_{res} (%)')

figure; hold on
histogram(X_res_per, "Normalization", "pdf");
histogram(X_res_per_sim, "Normalization", "pdf"); legend('experiment', 'simulation')
xlabel('X_{res} (%)'); ylabel('PDF')

%% Combustion Efficiency as function of AFR

% Combustion efficiency in percentage (%)
eta_c_per = eta_c' * 100;

% State
Mass_mg = [M_fuel, M_air] * 1e6; 

% Estimate the mean and covariance matrix
eta_c_mu    = mean([eta_c_per, Mass_mg]);
eta_c_Sigma = cov([eta_c_per,Mass_mg]);

% Conditional Gaussian
eta_c_per_sim = conditional_Gauss(eta_c_mu, eta_c_Sigma, Mass_mg);

figure; hold on
scatter(M_air./M_fuel, eta_c_per); ylim([0, 100]);
scatter(M_air./M_fuel, eta_c_per_sim); legend('experiment', 'simulation')
xlabel('Air-to-Fuel Ratio (-)'); ylabel('\eta_{c} (%)')

figure; hold on
histogram(eta_c_per, "Normalization", "pdf")
histogram(eta_c_per_sim, "Normalization", "pdf"); legend('experiment', 'simulation')
xlabel('\eta_{c} (%)'); ylabel('PDF')

%% Simluate entire file

% Initialize variables
M_fuel_sim = zeros(n_cycles,1);
M_air_sim = zeros(n_cycles,1);
eta_c_sim = zeros(n_cycles,1);
X_res_sim = zeros(n_cycles,1);
Q_gross_sim = zeros(n_cycles,1);

% Initial condition
M_fuel_sim(1) = M_fuel(1);
M_air_sim(1) = M_air(1);

for i = 1:n_cycles
    
    % State
    state = [M_fuel_sim(i); M_air_sim(i)];
    
    % Combustion efficiency
    eta_c_sim(i) = conditional_Gauss(eta_c_mu, eta_c_Sigma, state'*1e6) / 100;
    
    % Gross heat release
    Q_gross_sim(i) = eta_c_sim(i) * M_fuel_sim(i) * Q_LHV_mix(i);

    % Residual gas fraction
    X_res_sim(i) = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim(i)) / 100;

    % Residual mass matrix
    Matrix_res = X_res_sim(i) * [1 - eta_c_sim(i), 0; eta_c_sim(i), 1];

    % Fresh fuel and air
    input = [DI_quantity_c2c(i) + PI_quantity_c2c(i); Air_quantity_c2c(i)];

    if i < n_cycles
        % Calculate next cycle
        next_state = Matrix_A * state + input;
        M_fuel_sim(i+1) = next_state(1);
        M_air_sim(i+1) = next_state(2);
    end
end

figure; hold on
histogram(Q_gross, "Normalization", "pdf")
histogram(Q_gross_sim, "Normalization", "pdf"); legend('experiment', 'simulation')
xlabel('Q_{Gross} (J)'); ylabel('PDF')

figure; hold on
histogram(X_res, "Normalization", "pdf")
histogram(X_res_sim, "Normalization", "pdf"); legend('experiment', 'simulation')
xlabel('X_{res} (-)'); ylabel('PDF')

%% Combustion phasing

CA50 = Cylinder_1_Cycle_Data.CA50;
DI_timing = Cylinder_1_Cycle_Data.Injection_1_SOI;

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

% function [x1_sim, mu_cond, sigma_cond] = estimate_PDF(x1, x2, dist, dof)
% 
% % Combine data into a single matrix
% data = [x1, x2];
% 
% % Estimate the mean and covariance matrix
% mu = mean(data);
% Sigma = cov(data);
% 
% if dist == "Gauss"
%     [x1_sim, mu_cond, sigma_cond] = conditional_Gauss(mu, Sigma, x2);
% elseif dist == "t"
%     x1_sim = conditional_t(dof, mu, Sigma, x2);
% end
% end

% function x1_sim = conditional_t(dof, mu, Sigma, x2)
% 
% % Conditional t-distribution
% % x1 | x2 ~ t_nu(mu_cond, sigma_cond)
% 
% % Size of condiitonal variance
% p2 = size(x2,2);
% 
% % Conditional degrees of freedom
% nu_cond = dof + p2;
% 
% % Conditional mean
% mu_cond = mu(1) + Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * (x2 - mu(2:end))';
% 
% % Conditional Gaussian variance
% sigma_cond_Gauss = Sigma(1, 1) - Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * Sigma(2:end, 1);
% 
% % Squared Mahalanobis distance
% % d2 = (x2 - mu(2:end)) * Sigma(2:end, 2:end)^-1 * (x2 - mu(2:end))';
% d2 = dot(x2 - mu(2:end), (x2 - mu(2:end)) * Sigma(2:end, 2:end)^-1, 2);
% 
% % Conditional t-distribution standard deviation
% sigma_cond = sqrt((dof + d2) / nu_cond * sigma_cond_Gauss);
% 
% % Simulate Gaussian component
% y = normrnd(0, sigma_cond);
% 
% % Simulate chi-squared component
% u = chi2rnd(nu_cond, length(y), 1);
% 
% % Simulate residual gas fraction in percentage
% x1_sim = y ./ sqrt(u/nu_cond) + mu_cond';
% 
% end