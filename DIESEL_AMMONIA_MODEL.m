%%
clear all
close all
clc

% Load data
myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/';
load([myDir, '20240305/20240305_1200_DI6_6_SOI41_PFI8_5_038.mat'])

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
% Resize low=speed signals to match cycle-to-cycle values
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

% Exhaust mass per cucle (kg)
m_exhaust_c2c = m_diesel_c2c + m_NH3_c2c + m_air_c2c;

%% Extract exhaust gas information

% Molar mass (g/mol)
H2O_molar_mass = 18.01528;
CO2_molar_mass = 44.009;
O2_molar_mass = 31.999;
NH3_molar_mass = 17.031;
NO_molar_mass = 30.01;
N2_molar_mass = 28.0134;

% Mole fraction
H2O_mole_fraction = LowSpeed.FTIR_H2O_25 / 100;
CO2_mole_fraction = LowSpeed.FTIR_CO2_20 / 100;
O2_mole_fraction  = LowSpeed.Exhaust_O2 / 100;
NH3_mole_fraction = LowSpeed.FTIR_NH3_10000 / 1000000;
NO_mole_fraction  = LowSpeed.FTIR_NO_10000 / 1000000;
N2_mole_fraction  = 1 - (H2O_mole_fraction + CO2_mole_fraction + ...
                         O2_mole_fraction + NH3_mole_fraction + ...
                         NO_mole_fraction);

% Average molar mass of mixture
M_avg = H2O_mole_fraction * H2O_molar_mass + ...
        CO2_mole_fraction * CO2_molar_mass + ...
        O2_mole_fraction * O2_molar_mass + ...
        NH3_mole_fraction * NH3_molar_mass + ...
        NO_mole_fraction * NO_molar_mass + ...
        N2_mole_fraction * N2_molar_mass;

% Mass fraction of NH3
NH3_mass_fraction = NH3_mole_fraction * NH3_molar_mass ./ M_avg;
NH3_mass_fraction_c2c = interp1(1:n_samples, NH3_mass_fraction, c2c_index)';

% Unburned NH3 per cycle (kg)
m_NH3_unburned_c2c = m_exhaust_c2c .* NH3_mass_fraction;

%% Solve for unknowns

i = 154;
X_k = X_res(i);
m_diesel_k = m_diesel_c2c(i);
m_NH3_k = m_NH3_c2c(i);
m_air_k = m_air_c2c(i);
Q_gross_k = Q_gross(i);
m_NH3_unburned_k = m_NH3_unburned_c2c(i);

par = [X_k, m_diesel_k, m_NH3_k, m_air_k, Q_gross_k, m_NH3_unburned_k];

fun = @(x)dual_fuel_ss(x,par);
x = fmincon(fun, ones(1,5), [], [], [], [], [m_diesel_k, m_NH3_k, m_air_k, 0, 0], ones(1,5))

function y = dual_fuel_ss(x, par)
    
    % Unpack variables
    M_diesel = x(1);
    M_NH3 = x(2);
    M_air = x(3);
    eta_diesel = x(4);
    eta_NH3 = x(5);

    % Unpack parameters
    X_k = par(1);
    m_diesel_k = par(2);
    m_NH3_k = par(3);
    m_air_k = par(4);
    Q_gross_k = par(5);
    m_NH3_unburned_k = par(6);

    % Model parameters
    AFR_diesel = 14.5;
    AFR_NH3    = 6.04;
    Q_LHV_diesel = 44.1e6;
    Q_LHV_NH3    = 18.6e6;

    % Matrix
    A = [1 - X_k * (1 - eta_diesel),            0,               0;
                    0,               1 - X_k * (1 - eta_NH3),    0;
           AFR_diesel * eta_diesel,     AFR_NH3 * eta_NH3,    1 - X_k;
          Q_LHV_diesel * eta_diesel,   Q_LHV_NH3 * eta_NH3,      0;
                    0,                       1 - eta_NH3         0   ];
    
    % In-cylinder mass
    m = [M_diesel; M_NH3; M_air];

    % RHS
    b = [m_diesel_k; m_NH3_k; m_air_k; Q_gross_k; m_NH3_unburned_k];

    % Equation
    eqn = A * m - b;
    y = sum(eqn.^2);

end