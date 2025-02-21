%%
clear all
close all
clc

% Load data
myDir = "/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/20240305/";
files = dir(myDir+"*.mat");
n_files = length(files);

for j = 1:length(files)

% Load files
load([files(j).folder, '/', files(j).name])

% Get file name without extension
baseName = erase(files(j).name, ".mat");
baseName = baseName(15:end-10);

%% Estimate residual gas fraction

% Model parameters
EVC = -355;
EVO = 165;
gamma = 1.3;

% Crank angle axis
CA_deg = linspace(-359.8, 360, 3600);

% Number of cycles in file
n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);

% Cylinder pressure (Pa)
if length(Cylinder_1_Synch_Data.Cylinder_Pressure)/3600 == n_cycles
    Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure, [3600, n_cycles]) * 10 ^ -3;
else
    Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure(1:end-3600), [3600, n_cycles]) * 10 ^ -3;
end

% Valve timings
[~, i_evc] = min(abs(CA_deg - EVC));
[~, i_evo] = min(abs(CA_deg - EVO));

% Calculate the X_res (-)
X_res = zeros(n_cycles, 1);
for i = 1:n_cycles
    X_res(i) =  (Volume.Volume(i_evc) / Volume.Volume(i_evo)) * (Pcyl_CA(i_evc, i) / Pcyl_CA(i_evo, i)) ^ (1 / gamma);
end

%% Gross heat release (J)
Q_gross = Cylinder_1_Cycle_Data.Gross_Heat_Release';

%% Mass flows

% Low speed data sampling
n_samples = length(LowSpeed.Speed);
% Resize low-peed signals to match cycle-to-cycle values
c2c_index = linspace(1, n_samples, n_cycles);

% Diesel fuel mass per cycle (kg)
m_diesel = LowSpeed.Pilot_1_Injection_Quantity * 1e-6;
m_diesel_c2c = interp1(1:n_samples, m_diesel, c2c_index)';

% Ammonia fuel mass per cycle (kg)
m_NH3 = 2 * LowSpeed.ISB_Fuel_Flow_PFI ./ LowSpeed.Speed * 1e-3;
m_NH3_c2c = interp1(1:n_samples, m_NH3, c2c_index)';

% Air mass per cycle (kg)
m_air = 2 * LowSpeed.Mass_Air_Flow ./ LowSpeed.Speed * 1e-3;
m_air_c2c = interp1(1:n_samples, m_air, c2c_index)';

%% Combustion efficiency (-)
eta_c_mea = LowSpeed.Combustion_Efficiency' ./ 100;
eta_c_mea_c2c = interp1(1:n_samples, eta_c_mea, c2c_index)';

%% Solve for unknowns

% Initialize variables
M_diesel = zeros(n_cycles,1);
M_NH3 = zeros(n_cycles,1);
M_air = zeros(n_cycles,1);
eta_c = zeros(n_cycles,1);

for k = 1:n_cycles

    % Measured quantities
    par = [X_res(k), m_diesel_c2c(k), m_NH3_c2c(k), m_air_c2c(k), Q_gross(k)];
    
    % Root finding
    fun = @(x)dual_fuel(x,par);
    options = optimoptions('fsolve', 'Display', 'off');
    x_root = fsolve(fun, zeros(1,4), options);

    % Save for model
    M_diesel(k) = x_root(1);
    M_NH3(k) = x_root(2);
    M_air(k) = x_root(3);
    eta_c(k) = x_root(4);
end

%% Plot results

figure
subplot(3,1,1); plot(M_diesel*1e6); ylabel('In-cylinder diesel (mg)')
subplot(3,1,2); plot(M_NH3*1e6); ylabel('In-cylinder ammonia (mg)')
subplot(3,1,3); plot(M_air*1e6); ylabel('In-cylinder air (mg)')
xlabel('Cycles');
print(['Dual_Fuel_Plots/', baseName, 'Mass'], '-dpng', '-r300');

figure; hold on; box on;
plot(eta_c);
plot(eta_c_mea_c2c, 'LineWidth', 2); legend('Estimated', 'Measured');
ylabel('Combustion efficiency'); xlabel('Cycle')
print(['Dual_Fuel_Plots/', baseName, 'eta_c'], '-dpng', '-r300');

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

print(['Dual_Fuel_Plots/', baseName, 'X_res_model'], '-dpng', '-r300');

%% Combustion Efficiency as function of AFR

% Combustion efficiency in percentage (%)
eta_per = eta_c * 100;

% State
Mass_mg = [M_diesel, M_NH3, M_air] * 1e6; 

% Estimate the mean and covariance matrix
eta_c_mu    = mean([eta_per, Mass_mg]);
eta_c_Sigma = cov([eta_per, Mass_mg]);

% Conditional Gaussian
eta_c_per_model = conditional_Gauss(eta_c_mu, eta_c_Sigma, Mass_mg);

figure
subplot(2,1,1); hold on; box on
scatter(M_air./(M_diesel + M_NH3), eta_per); ylim([0, 100]);
scatter(M_air./(M_diesel + M_NH3), eta_c_per_model); legend('Experimental', 'Conditional Gaussian')
xlabel('Air-to-Fuel Ratio (-)'); ylabel('Combustion Efficiency (%)')

subplot(2,1,2); hold on; box on
histogram(eta_per, "Normalization", "pdf")
histogram(eta_c_per_model, "Normalization", "pdf");
legend('Experimental', 'Conditional Gaussian')
xlabel('Combustion Efficiency (%)'); ylabel('PDF')

print(['Dual_Fuel_Plots/', baseName, 'eta_c_model'], '-dpng', '-r300');

%% Combustion phasing
CA50 = Cylinder_1_Cycle_Data.CA50';
DI_timing = -Cylinder_1_Cycle_Data.Injection_1_SOI';

%% Desired operating conditions

DI_quantity_des = mean(LowSpeed.Pilot_1_Injection_Quantity);
DI_timing_des = -mean(Cylinder_1_Cycle_Data.Injection_1_SOI);

%% CA50 as function of Injection timing

% State
state_ca50 = [M_diesel* 1e6, M_NH3* 1e6, M_air*1e4, DI_timing] ; 

% Estimate the mean and covariance matrix
CA50_mu    = mean([CA50, state_ca50]);
CA50_Sigma = cov([CA50, state_ca50]);

% Conditional Gaussian
CA50_model = conditional_Gauss(CA50_mu, CA50_Sigma, state_ca50);

figure
subplot(2,1,1); hold on; box on
scatter(M_air./(M_diesel + M_NH3), CA50);
scatter(M_air./(M_diesel + M_NH3), CA50_model); legend('Experimental', 'Conditional Gaussian')
xlabel('Air-to-Fuel Ratio (-)'); ylabel('CA50 (deg)')

subplot(2,1,2); hold on; box on
histogram(CA50, "Normalization", "pdf")
histogram(CA50_model, "Normalization", "pdf");
legend('Experimental', 'Conditional Gaussian')
xlabel('CA50 (deg)'); ylabel('PDF')

print(['Dual_Fuel_Plots/', baseName, 'CA50_model'], '-dpng', '-r300');

%% Simulate dynamic system

% Initialize variables
M_diesel_sim = zeros(n_cycles,1);
M_NH3_sim = zeros(n_cycles,1);
M_air_sim = zeros(n_cycles,1);
eta_c_sim = zeros(n_cycles,1);
X_res_sim = zeros(n_cycles,1);
Q_gross_sim = zeros(n_cycles,1);
CA50_sim = zeros(n_cycles,1);

% Initial condition
M_diesel_sim(1) = M_diesel(1);
M_NH3_sim(1) = M_NH3(1);
M_air_sim(1) = M_air(1);

% Inputs
Diesel_fuel = m_diesel_c2c;
NH3_fuel = m_NH3_c2c;
Fresh_air = m_air_c2c;
Injection_time = DI_timing;

% Model parameters
AFR_diesel = 14.5;
AFR_NH3    = 6.04;
Q_LHV_diesel = 44.1e6;
Q_LHV_NH3    = 18.6e6;

for i = 1:n_cycles

    % State
    state = [M_diesel_sim(i); M_NH3_sim(i); M_air_sim(i)];
    state_CA50 = [M_diesel_sim(i) * 1e6, M_NH3_sim(i) * 1e6, ...
                  M_air_sim(i) * 1e4, Injection_time(i)];

    % CA50
    CA50_sim(i) = conditional_Gauss(CA50_mu, CA50_Sigma, state_CA50);
    
    % Combustion efficiency
    eta_c_sim(i) = conditional_Gauss(eta_c_mu, eta_c_Sigma, state'*1e6) / 100;

    % Gross heat release
    Q_gross_sim(i) = eta_c_sim(i) * (Q_LHV_diesel * M_diesel_sim(i) + Q_LHV_NH3 * M_NH3_sim(i));

    % Residual gas fraction
    X_res_sim(i) = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim(i)) / 100;

    % Residual mass matrix
    Matrix_res = X_res_sim(i) * [     1 - eta_c_sim(i),                0,           0
                                              0,               1 -  eta_c_sim(i),   0
                                -AFR_diesel * eta_c_sim(i), -AFR_NH3 * eta_c_sim(i), 1];

    % Fresh fuel and air
    input = [Diesel_fuel(i); NH3_fuel(i); Fresh_air(i)];

    if i < n_cycles
        % Calculate next cycle
        next_state = Matrix_res * state + input;
        M_diesel_sim(i+1) = next_state(1);
        M_NH3_sim(i+1) = next_state(2);
        M_air_sim(i+1) = next_state(3);
    end
end

%% Fix initial condition
M_diesel_sim(1) = M_diesel_sim(2);
M_NH3_sim(1) = M_NH3_sim(2);
M_air_sim(1) = M_air_sim(2);

figure
subplot(6,2,1); hold on; box on;
plot(M_diesel*1e6); plot(M_diesel_sim*1e6); ylabel('M_{diesel} (mg)');
xLimits = ylim; title('Timeseries')
subplot(6,2,2); hold on; box on; xlim(xLimits)
histogram(M_diesel*1e6, "Normalization", "pdf")
histogram(M_diesel_sim*1e6, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse'); title('PDF')

subplot(6,2,3); hold on; box on;
plot(M_NH3*1e6); plot(M_NH3_sim*1e6); ylabel('M_{NH3} (mg)');
xLimits = ylim;
subplot(6,2,4); hold on; box on; xlim(xLimits)
histogram(M_NH3*1e6, "Normalization", "pdf")
histogram(M_NH3_sim*1e6, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

subplot(6,2,5); hold on; box on;
plot(M_air*1e6); plot(M_air_sim*1e6); ylabel('M_{air} (mg)'); xLimits = ylim;
subplot(6,2,6); hold on; box on; xlim(xLimits)
histogram(M_air*1e6, "Normalization", "pdf")
histogram(M_air_sim*1e6, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

subplot(6,2,7); hold on; box on;
plot(Q_gross); plot(Q_gross_sim); ylabel('Q_{gross} (J)'); xLimits = ylim;
subplot(6,2,8); hold on; box on; xlim(xLimits)
histogram(Q_gross, "Normalization", "pdf")
histogram(Q_gross_sim, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

subplot(6,2,9); hold on; box on;
plot(X_res*100); plot(X_res_sim*100); ylabel('X_{res} (%)'); xLimits = ylim;
subplot(6,2,10); hold on; box on; xlim(xLimits)
histogram(X_res*100, "Normalization", "pdf")
histogram(X_res_sim*100, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

subplot(6,2,11); hold on; box on;
plot(CA50); plot(CA50_sim); ylabel('CA50 (deg)'); xLimits = ylim;
subplot(6,2,12); hold on; box on; xlim(xLimits)
histogram(CA50, "Normalization", "pdf")
histogram(CA50_sim, "Normalization", "pdf");
legend('Data', 'Simulator','Location','eastoutside')
view(90, 90); set(gca, 'XDir', 'reverse');

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.7]);
print(['Dual_Fuel_Plots/', baseName, 'Simulator'], '-dpng', '-r300');

%% Safe file

model_data = table(M_diesel_sim, M_NH3_sim, M_air_sim, ...
                   Diesel_fuel,  NH3_fuel,  Fresh_air,...
                   Q_gross, CA50, DI_timing, ...
                   X_res);
writetable(model_data, ['Dual_Fuel_Data/', baseName, 'Data.csv']);
save(['Dual_Fuel_Data/', baseName, 'Parameters.mat'], ...
    'eta_c_mu', 'eta_c_Sigma', ...
    'X_res_mu', 'X_res_Sigma', ...
    'CA50_mu', 'CA50_Sigma', ...
    'Q_LHV_diesel', 'Q_LHV_NH3', ...
    'DI_quantity_des', 'DI_timing_des');

end

%% Functions
function y = dual_fuel(x, par)
    
    % Unpack variables
    M_diesel = x(1);
    M_NH3 = x(2);
    M_air = x(3);
    eta = x(4);

    % Unpack parameters
    X_k = par(1);
    m_diesel_k = par(2);
    m_NH3_k = par(3);
    m_air_k = par(4);
    Q_gross_k = par(5);

    % Model parameters
    AFR_diesel = 14.5;
    AFR_NH3    = 6.04;
    Q_LHV_diesel = 44.1e6;
    Q_LHV_NH3    = 18.6e6;

    % Matrix
    A = [ 1 - X_k * (1 - eta),            0,            0;
                    0,           1 - X_k * (1 - eta),   0;
         X_k * AFR_diesel * eta, X_k * AFR_NH3 * eta, 1 - X_k;
           Q_LHV_diesel * eta,     Q_LHV_NH3 * eta,     0];
    
    % In-cylinder mass
    m = [M_diesel; M_NH3; M_air];

    % RHS
    b = [m_diesel_k; m_NH3_k; m_air_k; Q_gross_k];

    % Equation
    y = A * m - b;
    y = y .* [1e6, 1e5, 1e3, 1e-2]';

end

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