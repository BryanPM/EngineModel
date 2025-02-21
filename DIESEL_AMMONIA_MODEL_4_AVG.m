%%
clear all
close all
clc

% Load data
myDir = "/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/20240305/";
files = dir(myDir+"*.mat");
n_files = length(files);

MEASUREMENTS = zeros(4,n_files);
ROOT_FIND = zeros(4,n_files);

for j = 1:length(files)
load([files(j).folder, '/', files(j).name])

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

%% Combustion efficiency (-)
eta_c = LowSpeed.Combustion_Efficiency' ./ 100;

%% Mass flows

% Low speed data sampling
n_samples = length(LowSpeed.Speed);
% Resize low-speed signals to match cycle-to-cycle values
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

%% Solve for unknowns

% Measured quantities
X_mean = mean(X_res);
m_diesel_mean = mean(m_diesel_c2c);
m_NH3_mean = mean(m_NH3_c2c);
m_air_mean = mean(m_air_c2c);
Q_gross_mean = mean(Q_gross);
par = [X_mean, m_diesel_mean, m_NH3_mean, m_air_mean, Q_gross_mean];
eta_c_mean = mean(eta_c);

% Root finding
fun = @(x)dual_fuel(x,par);
options = optimoptions('fsolve', 'Display', 'off');
x_root = fsolve(fun, zeros(1,4), options);

% Print
FILE_NAME{j} = files(j).name;
MEASUREMENTS(:,j) = [m_diesel_mean * 1e6;
                     m_NH3_mean * 1e6;
                     m_air_mean * 1e3;
                     eta_c_mean * 1e2];
ROOT_FIND(:,j) = x_root;

end

%% Plots

figure(1); clf; hold on
plot(MEASUREMENTS(1,:), 'o', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'r')
plot(ROOT_FIND(1,:) * 1e6, 's', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'g')
legend({'m_{in}', 'M_{fuel}'}, 'Location', 'best');
xlabel('File number')
ylabel('Diesel fuel (mg)')

figure(2); clf; hold on
plot(MEASUREMENTS(2,:), 'o', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'r')
plot(ROOT_FIND(2,:) * 1e6, 's', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'g')
legend({'m_{in}', 'M_{fuel}'}, 'Location', 'best');
xlabel('File number')
ylabel('Ammonia fuel (mg)')

figure(3); clf; hold on
plot(MEASUREMENTS(3,:), 'o', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'r')
plot(ROOT_FIND(3,:) * 1e3, 's', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'g')
legend({'m_{in}', 'M_{air}'}, 'Location', 'best');
xlabel('File number')
ylabel('Air mass (g)')

figure(4); clf; hold on
plot(MEASUREMENTS(4,:), 'o', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'r')
plot(ROOT_FIND(4,:) * 1e2, 's', 'LineWidth', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'g')
legend({'Measured', 'Estimated'}, 'Location', 'best');
xlabel('File number')
ylabel('Combustion Efficiency (-)')

%% Functions
function y = dual_fuel(x, par)
    
    % Unpack variables
    M_diesel = x(1);
    M_NH3 = x(2);
    M_air = x(3);
    eta_c = x(4);

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
    A = [ 1 - X_k * (1 - eta_c),               0,            0;
                    0,             1 - X_k * (1 - eta_c),    0;
         X_k * AFR_diesel * eta_c, X_k * AFR_NH3 * eta_c, 1 - X_k;
           Q_LHV_diesel * eta_c,     Q_LHV_NH3 * eta_c,      0];
    
    % In-cylinder mass
    m = [M_diesel; M_NH3; M_air];

    % RHS
    b = [m_diesel_k; m_NH3_k; m_air_k; Q_gross_k];

    % Equation
    y = A * m - b;
    y = y .* [1e6, 1e5, 1e3, 1e-2]';

end
