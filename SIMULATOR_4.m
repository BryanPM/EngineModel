%% Tabula Rasa
clear all
close all
clc

%% Measured data

% Load data
myDir = 'Dual_Fuel_Data/';
myFiles = dir(fullfile(myDir, 'DI*SOI*_Data.csv'));

% Initialize variables
Q_gross_all = [];
m_diesel_all = [];
m_NH3_all = [];
m_air_all = [];
CA50_all = [];
DI_SOI_all = [];

% Cycles per file
n_cycles_file = 100;

for i = 1:length(myFiles)

% Load csv files
fileName = myFiles(i).name;
cycleData = readmatrix([myDir, fileName]);

% Extract initial condition
if i == 1
    M_diesel_init = cycleData(1,1);
    M_NH3_init = cycleData(1,2);
    M_air_init  = cycleData(1,3);
end

% Extract time series
m_diesel_all = [m_diesel_all; cycleData(1:n_cycles_file,4)];
m_NH3_all = [m_NH3_all; cycleData(1:n_cycles_file,5)];
m_air_all = [m_air_all; cycleData(1:n_cycles_file,6)];
Q_gross_all = [Q_gross_all; cycleData(1:n_cycles_file,7)];
CA50_all = [CA50_all; cycleData(1:n_cycles_file,8)];
DI_SOI_all = [DI_SOI_all; cycleData(1:n_cycles_file,9)];

end

%% Simulator

% Total number of cycles
n_cycles = length(Q_gross_all);

% Initialize variables
M_diesel_sim = zeros(n_cycles,1);
M_NH3_sim = zeros(n_cycles,1);
M_air_sim = zeros(n_cycles,1);
eta_c_sim = zeros(n_cycles,1);
X_res_sim = zeros(n_cycles,1);
Q_gross_sim = zeros(n_cycles,1);
CA50_sim = zeros(n_cycles,1);

% Initial condition
M_diesel_sim(1) = M_diesel_init;
M_NH3_sim(1) = M_NH3_init;
M_air_sim(1)  = M_air_init;

% Model parameters
AFR_diesel = 14.5;
AFR_NH3    = 6.04;
Q_LHV_diesel = 44.1e6;
Q_LHV_NH3    = 18.6e6;

% Import lookup tables
lookup_tables = load_lookup_tables(myDir);

for i = 1:n_cycles

    % Mass inputs
    Diesel_fuel = m_diesel_all(i);
    NH3_fuel = m_NH3_all(i);
    Fresh_air = m_air_all(i);

    % Other inputs
    Diesel_SOI = DI_SOI_all(i);

    
    % State
    state = [M_diesel_sim(i); M_NH3_sim(i); M_air_sim(i)];
    state_CA50 = [M_diesel_sim(i) * 1e6, M_NH3_sim(i) * 1e6, ...
                  M_air_sim(i) * 1e4, Diesel_SOI];
    
    % Conditional distribution parameters
    [eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, CA50_mu, ...
    CA50_Sigma] = Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables);

    % Simulate CA50
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
    input = [Diesel_fuel; NH3_fuel; Fresh_air];

    if i < n_cycles
        % Calculate next cycle
        next_state = Matrix_res * state + input;
        M_diesel_sim(i+1) = next_state(1);
        M_NH3_sim(i+1) = next_state(2);
        M_air_sim(i+1) = next_state(3);
    end
end

% Plot
RMSE = sqrt(mean((Q_gross_all - Q_gross_sim).^2));
figure; clf; hold on; box on;
plot(Q_gross_all)
plot(Q_gross_sim)
ylabel('Q_{gross} (J)'); legend('Experimental', 'Simulated')
title("RMSE = " + num2str(RMSE, 3) + " (J)")
screenSize = get(0, 'ScreenSize'); xlabel('Cycles');
originalHeight = 400; % You can adjust this height value as needed
set(gcf, 'Units', 'pixels', 'Position', [0, 0, screenSize(3), originalHeight]);

RMSE = sqrt(mean((CA50_all - CA50_sim).^2));
figure; clf; hold on; box on;
plot(CA50_all)
plot(CA50_sim)
ylabel('CA50 (aTDC)'); legend('Experimental', 'Simulated')
title("RMSE = " + num2str(RMSE, 2) + " (deg)")
screenSize = get(0, 'ScreenSize'); xlabel('Cycles');
originalHeight = 400; % You can adjust this height value as needed
set(gcf, 'Units', 'pixels', 'Position', [0, 0, screenSize(3), originalHeight]);

%% Additional functions
function x1_sim = conditional_Gauss(mu, Sigma, x2)

% Conditional Gaussian
% x1 | x2 ~ N(mu_cond, sigma_cond)

% Conditional standard deviation
sigma_cond = sqrt(Sigma(1, 1) - Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * Sigma(2:end, 1));
if sigma_cond < 0
    sigma_cond = 0;
    warning('Negative conditional variance, setting it to zero')
end
% Gaussian Conditional mean
mu_cond = mu(1) + Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * (x2 - mu(2:end))';

% Simulate residual gas fraction in percentage
x1_sim = normrnd(mu_cond, sigma_cond);

end

function [eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, CA50_mu, ...
    CA50_Sigma] = Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables)

% Joint mean for combustion efficiency distribution
mu_eta_1_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_eta_1, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_eta_2_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_eta_2, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_eta_3_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_eta_3, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_eta_4_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_eta_4, ...
                        Diesel_fuel*1e6, Diesel_SOI);

% Covariance matrix for combustion efficiency distribution
Sigma_eta_11_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_11, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_12_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_12, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_13_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_13, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_14_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_14, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_22_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_22, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_23_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_23, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_24_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_24, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_33_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_33, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_34_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_34, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_44_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_44, ...
                            Diesel_fuel*1e6, Diesel_SOI);

% Joint mean for residual gas fraction
mu_X_1_eval = interp2(lookup_tables.DI_QTY_interp, ...
                      lookup_tables.DI_SOI_interp, ...
                      lookup_tables.mu_X_1, ...
                      Diesel_fuel*1e6, Diesel_SOI);
mu_X_2_eval = interp2(lookup_tables.DI_QTY_interp, ...
                      lookup_tables.DI_SOI_interp, ...
                      lookup_tables.mu_X_2, ...
                      Diesel_fuel*1e6, Diesel_SOI);

% Covariance matrix for residual gas fraction
Sigma_X_11_eval = interp2(lookup_tables.DI_QTY_interp, ...
                          lookup_tables.DI_SOI_interp, ...
                          lookup_tables.Sigma_X_11, ...
                          Diesel_fuel*1e6, Diesel_SOI);
Sigma_X_12_eval = interp2(lookup_tables.DI_QTY_interp, ...
                          lookup_tables.DI_SOI_interp, ...
                          lookup_tables.Sigma_X_12, ...
                          Diesel_fuel*1e6, Diesel_SOI);
Sigma_X_22_eval = interp2(lookup_tables.DI_QTY_interp, ...
                          lookup_tables.DI_SOI_interp, ...
                          lookup_tables.Sigma_X_22, ...
                          Diesel_fuel*1e6, Diesel_SOI);

% Joint mean for CA50 distribution
mu_CA50_1_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_CA50_1, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_CA50_2_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_CA50_2, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_CA50_3_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_CA50_3, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_CA50_4_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_CA50_4, ...
                        Diesel_fuel*1e6, Diesel_SOI);
mu_CA50_5_eval = interp2(lookup_tables.DI_QTY_interp, ...
                        lookup_tables.DI_SOI_interp, ...
                        lookup_tables.mu_CA50_5, ...
                        Diesel_fuel*1e6, Diesel_SOI);

% Covariance matrix for CA50 distribution
Sigma_CA50_11_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_11, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_12_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_12, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_13_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_13, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_14_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_14, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_15_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_15, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_22_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_22, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_23_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_23, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_24_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_24, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_25_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_25, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_33_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_33, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_34_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_34, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_35_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_35, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_44_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_44, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_45_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_45, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_CA50_55_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_CA50_55, ...
                            Diesel_fuel*1e6, Diesel_SOI);

% Return values
X_res_mu    = [mu_X_1_eval, mu_X_2_eval];
X_res_Sigma = [Sigma_X_11_eval, Sigma_X_12_eval;
               Sigma_X_12_eval, Sigma_X_22_eval];

% Return values
eta_c_mu    = [mu_eta_1_eval, mu_eta_2_eval, mu_eta_3_eval, mu_eta_4_eval];
eta_c_Sigma = [Sigma_eta_11_eval, Sigma_eta_12_eval, Sigma_eta_13_eval, Sigma_eta_14_eval;
               Sigma_eta_12_eval, Sigma_eta_22_eval, Sigma_eta_23_eval, Sigma_eta_24_eval;
               Sigma_eta_13_eval, Sigma_eta_23_eval, Sigma_eta_33_eval, Sigma_eta_34_eval;
               Sigma_eta_14_eval, Sigma_eta_24_eval, Sigma_eta_34_eval, Sigma_eta_44_eval];

% Return values
CA50_mu    = [mu_CA50_1_eval, mu_CA50_2_eval, mu_CA50_3_eval, mu_CA50_4_eval, mu_CA50_5_eval];
CA50_Sigma = [Sigma_CA50_11_eval, Sigma_CA50_12_eval, Sigma_CA50_13_eval, Sigma_CA50_14_eval, Sigma_CA50_15_eval;
              Sigma_CA50_12_eval, Sigma_CA50_22_eval, Sigma_CA50_23_eval, Sigma_CA50_24_eval, Sigma_CA50_25_eval;
              Sigma_CA50_13_eval, Sigma_CA50_23_eval, Sigma_CA50_33_eval, Sigma_CA50_34_eval, Sigma_CA50_35_eval;
              Sigma_CA50_14_eval, Sigma_CA50_24_eval, Sigma_CA50_34_eval, Sigma_CA50_44_eval, Sigma_CA50_45_eval;
              Sigma_CA50_15_eval, Sigma_CA50_25_eval, Sigma_CA50_35_eval, Sigma_CA50_45_eval, Sigma_CA50_55_eval];

end

function lookup_tables = load_lookup_tables(myDir)

% Import data
mu_eta_1 = readmatrix([myDir, 'mu_eta_1']);
mu_eta_2 = readmatrix([myDir, 'mu_eta_2']);
mu_eta_3 = readmatrix([myDir, 'mu_eta_3']);
mu_eta_4 = readmatrix([myDir, 'mu_eta_4']);
Sigma_eta_11 = readmatrix([myDir, 'Sigma_eta_11']);
Sigma_eta_12 = readmatrix([myDir, 'Sigma_eta_12']);
Sigma_eta_13 = readmatrix([myDir, 'Sigma_eta_13']);
Sigma_eta_14 = readmatrix([myDir, 'Sigma_eta_14']);
Sigma_eta_22 = readmatrix([myDir, 'Sigma_eta_22']);
Sigma_eta_23 = readmatrix([myDir, 'Sigma_eta_23']);
Sigma_eta_24 = readmatrix([myDir, 'Sigma_eta_24']);
Sigma_eta_33 = readmatrix([myDir, 'Sigma_eta_33']);
Sigma_eta_34 = readmatrix([myDir, 'Sigma_eta_34']);
Sigma_eta_44 = readmatrix([myDir, 'Sigma_eta_44']);
mu_X_1 = readmatrix([myDir, 'mu_X_1']);
mu_X_2 = readmatrix([myDir, 'mu_X_2']);
Sigma_X_11 = readmatrix([myDir, 'Sigma_X_11']);
Sigma_X_12 = readmatrix([myDir, 'Sigma_X_12']);
Sigma_X_22 = readmatrix([myDir, 'Sigma_X_22']);
mu_CA50_1 = readmatrix([myDir, 'mu_CA50_1']);
mu_CA50_2 = readmatrix([myDir, 'mu_CA50_2']);
mu_CA50_3 = readmatrix([myDir, 'mu_CA50_3']);
mu_CA50_4 = readmatrix([myDir, 'mu_CA50_4']);
mu_CA50_5 = readmatrix([myDir, 'mu_CA50_5']);
Sigma_CA50_11 = readmatrix([myDir, 'Sigma_CA50_11']);
Sigma_CA50_12 = readmatrix([myDir, 'Sigma_CA50_12']);
Sigma_CA50_13 = readmatrix([myDir, 'Sigma_CA50_13']);
Sigma_CA50_14 = readmatrix([myDir, 'Sigma_CA50_14']);
Sigma_CA50_15 = readmatrix([myDir, 'Sigma_CA50_15']);
Sigma_CA50_22 = readmatrix([myDir, 'Sigma_CA50_22']);
Sigma_CA50_23 = readmatrix([myDir, 'Sigma_CA50_23']);
Sigma_CA50_24 = readmatrix([myDir, 'Sigma_CA50_24']);
Sigma_CA50_25 = readmatrix([myDir, 'Sigma_CA50_25']);
Sigma_CA50_33 = readmatrix([myDir, 'Sigma_CA50_33']);
Sigma_CA50_34 = readmatrix([myDir, 'Sigma_CA50_34']);
Sigma_CA50_35 = readmatrix([myDir, 'Sigma_CA50_35']);
Sigma_CA50_44 = readmatrix([myDir, 'Sigma_CA50_44']);
Sigma_CA50_45 = readmatrix([myDir, 'Sigma_CA50_45']);
Sigma_CA50_55 = readmatrix([myDir, 'Sigma_CA50_55']);
DI_QTY_interp = readmatrix([myDir, 'DI_QTY_interp']);
DI_SOI_interp = readmatrix([myDir, 'DI_SOI_interp']);

% Create datebase
lookup_tables.mu_eta_1 = mu_eta_1;
lookup_tables.mu_eta_2 = mu_eta_2;
lookup_tables.mu_eta_3 = mu_eta_3;
lookup_tables.mu_eta_4 = mu_eta_4;
lookup_tables.Sigma_eta_11 = Sigma_eta_11;
lookup_tables.Sigma_eta_12 = Sigma_eta_12;
lookup_tables.Sigma_eta_13 = Sigma_eta_13;
lookup_tables.Sigma_eta_14 = Sigma_eta_14;
lookup_tables.Sigma_eta_22 = Sigma_eta_22;
lookup_tables.Sigma_eta_23 = Sigma_eta_23;
lookup_tables.Sigma_eta_24 = Sigma_eta_24;
lookup_tables.Sigma_eta_33 = Sigma_eta_33;
lookup_tables.Sigma_eta_34 = Sigma_eta_34;
lookup_tables.Sigma_eta_44 = Sigma_eta_44;
lookup_tables.mu_X_1 = mu_X_1;
lookup_tables.mu_X_2 = mu_X_2;
lookup_tables.Sigma_X_11 = Sigma_X_11;
lookup_tables.Sigma_X_12 = Sigma_X_12;
lookup_tables.Sigma_X_22 = Sigma_X_22;
lookup_tables.mu_CA50_1 = mu_CA50_1;
lookup_tables.mu_CA50_2 = mu_CA50_2;
lookup_tables.mu_CA50_3 = mu_CA50_3;
lookup_tables.mu_CA50_4 = mu_CA50_4;
lookup_tables.mu_CA50_5 = mu_CA50_5;
lookup_tables.Sigma_CA50_11 = Sigma_CA50_11;
lookup_tables.Sigma_CA50_12 = Sigma_CA50_12;
lookup_tables.Sigma_CA50_13 = Sigma_CA50_13;
lookup_tables.Sigma_CA50_14 = Sigma_CA50_14;
lookup_tables.Sigma_CA50_15 = Sigma_CA50_15;
lookup_tables.Sigma_CA50_22 = Sigma_CA50_22;
lookup_tables.Sigma_CA50_23 = Sigma_CA50_23;
lookup_tables.Sigma_CA50_24 = Sigma_CA50_24;
lookup_tables.Sigma_CA50_25 = Sigma_CA50_25;
lookup_tables.Sigma_CA50_33 = Sigma_CA50_33;
lookup_tables.Sigma_CA50_34 = Sigma_CA50_34;
lookup_tables.Sigma_CA50_35 = Sigma_CA50_35;
lookup_tables.Sigma_CA50_44 = Sigma_CA50_44;
lookup_tables.Sigma_CA50_45 = Sigma_CA50_45;
lookup_tables.Sigma_CA50_55 = Sigma_CA50_55;
lookup_tables.DI_QTY_interp = DI_QTY_interp;
lookup_tables.DI_SOI_interp = DI_SOI_interp;

end