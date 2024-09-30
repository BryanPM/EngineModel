%% Tabula Rasa
clear all
close all
clc

%% Measured data

% Load data
myDir = 'Model_Data/';
myFiles = dir(fullfile(myDir, 'DI*SOI*_Data.csv'));

% Initialize variables
Q_gross_all = [];
m_diesel_all = [];
m_ammonia_all = [];
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
    M_fuel_init = cycleData(1,1);
    M_air_init  = cycleData(2,2);
end

% Extract time series
Q_gross_all = [Q_gross_all; cycleData(1:n_cycles_file,4)];
m_diesel_all = [m_diesel_all; cycleData(1:n_cycles_file,6)];
m_ammonia_all = [m_ammonia_all; cycleData(1:n_cycles_file,7)];
m_air_all = [m_air_all; cycleData(1:n_cycles_file,8)];
CA50_all = [CA50_all; cycleData(1:n_cycles_file,9)];
DI_SOI_all = [DI_SOI_all; cycleData(1:n_cycles_file,10)];

end

%% Simulator

% Lower heating values (J)
Q_LHV_diesel  = 44.1 * 1e6;
Q_LHV_ammonia = 18.6 * 1e6;

% Total number of cycles
n_cycles = length(Q_gross_all);

% Initialize variables
M_fuel_sim = zeros(n_cycles,1);
M_air_sim = zeros(n_cycles,1);
eta_c_sim = zeros(n_cycles,1);
X_res_sim = zeros(n_cycles,1);
Q_gross_sim = zeros(n_cycles,1);
CA50_sim = zeros(n_cycles,1);

% Initial condition
M_fuel_sim(1) = M_fuel_init;
M_air_sim(1)  = M_air_init;

% Import lookup tables
lookup_tables = load_lookup_tables(myDir);

for i = 1:n_cycles

    % Diesel inputs
    Diesel_SOI  = DI_SOI_all(i);
    Diesel_fuel = m_diesel_all(i);

    % Other inputs
    Ammonia_fuel = m_ammonia_all(i);
    Fresh_air = m_air_all(i);
    
    % State
    state = [M_fuel_sim(i); M_air_sim(i)];
    
    % Conditional distribution parameters
    [eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, mu_CA50_eval, ...
    Sigma_CA50_eval] = Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables);

    % Simulate CA50
    CA50_sim(i) = normrnd(mu_CA50_eval, Sigma_CA50_eval);

    % Combustion efficiency
    eta_c_sim(i) = conditional_Gauss(eta_c_mu, eta_c_Sigma, state'*1e6) / 100;
    
    % Effective LHV
    Q_LHV_eff = (Diesel_fuel * Q_LHV_diesel + Ammonia_fuel * Q_LHV_ammonia) ./ (Diesel_fuel + Ammonia_fuel);

    % Gross heat release
    Q_gross_sim(i) = eta_c_sim(i) * M_fuel_sim(i) * Q_LHV_eff;

    % Residual gas fraction
    X_res_sim(i) = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim(i)) / 100;

    % Residual mass matrix
    Matrix_res = X_res_sim(i) * [1 - eta_c_sim(i), 0; eta_c_sim(i), 1];

    % Fresh fuel and air
    input = [Diesel_fuel + Ammonia_fuel; Fresh_air];

    if i < n_cycles
        % Calculate next cycle
        next_state = Matrix_res * state + input;
        M_fuel_sim(i+1) = next_state(1);
        M_air_sim(i+1) = next_state(2);
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

function [eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, mu_CA50_eval, ...
          Sigma_CA50_eval] = Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables)

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
Sigma_eta_22_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_22, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_23_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_23, ...
                            Diesel_fuel*1e6, Diesel_SOI);
Sigma_eta_33_eval = interp2(lookup_tables.DI_QTY_interp, ...
                            lookup_tables.DI_SOI_interp, ...
                            lookup_tables.Sigma_eta_33, ...
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

% Mean for CA50 distribution
mu_CA50_eval = interp2(lookup_tables.DI_QTY_interp, ...
                       lookup_tables.DI_SOI_interp, ...
                       lookup_tables.mu_CA50, ...
                       Diesel_fuel*1e6, Diesel_SOI);

% Covariance matrix for CA50 distribution
Sigma_CA50_eval = interp2(lookup_tables.DI_QTY_interp, ...
                          lookup_tables.DI_SOI_interp, ...
                          lookup_tables.Sigma_CA50, ...
                          Diesel_fuel*1e6, Diesel_SOI);

% Return values
X_res_mu    = [mu_X_1_eval, mu_X_2_eval];
X_res_Sigma = [Sigma_X_11_eval, Sigma_X_12_eval;
               Sigma_X_12_eval, Sigma_X_22_eval];

% Return values
eta_c_mu    = [mu_eta_1_eval, mu_eta_2_eval, mu_eta_3_eval];
eta_c_Sigma = [Sigma_eta_11_eval, Sigma_eta_12_eval, Sigma_eta_13_eval;
               Sigma_eta_12_eval, Sigma_eta_22_eval, Sigma_eta_23_eval;
               Sigma_eta_13_eval, Sigma_eta_23_eval, Sigma_eta_33_eval];

end

function lookup_tables = load_lookup_tables(myDir)

% Import data
mu_eta_1 = readmatrix([myDir, 'mu_eta_1']);
mu_eta_2 = readmatrix([myDir, 'mu_eta_2']);
mu_eta_3 = readmatrix([myDir, 'mu_eta_3']);
Sigma_eta_11 = readmatrix([myDir, 'Sigma_eta_11']);
Sigma_eta_12 = readmatrix([myDir, 'Sigma_eta_12']);
Sigma_eta_13 = readmatrix([myDir, 'Sigma_eta_13']);
Sigma_eta_22 = readmatrix([myDir, 'Sigma_eta_22']);
Sigma_eta_23 = readmatrix([myDir, 'Sigma_eta_23']);
Sigma_eta_33 = readmatrix([myDir, 'Sigma_eta_33']);
mu_X_1 = readmatrix([myDir, 'mu_X_1']);
mu_X_2 = readmatrix([myDir, 'mu_X_2']);
Sigma_X_11 = readmatrix([myDir, 'Sigma_X_11']);
Sigma_X_12 = readmatrix([myDir, 'Sigma_X_12']);
Sigma_X_22 = readmatrix([myDir, 'Sigma_X_22']);
mu_CA50 = readmatrix([myDir, 'mu_CA50']);
Sigma_CA50 = readmatrix([myDir, 'Sigma_CA50']);
DI_QTY_interp = readmatrix([myDir, 'DI_QTY_interp']);
DI_SOI_interp = readmatrix([myDir, 'DI_SOI_interp']);

% Create datebase
lookup_tables.mu_eta_1 = mu_eta_1;
lookup_tables.mu_eta_2 = mu_eta_2;
lookup_tables.mu_eta_3 = mu_eta_3;
lookup_tables.Sigma_eta_11 = Sigma_eta_11;
lookup_tables.Sigma_eta_12 = Sigma_eta_12;
lookup_tables.Sigma_eta_13 = Sigma_eta_13;
lookup_tables.Sigma_eta_22 = Sigma_eta_22;
lookup_tables.Sigma_eta_23 = Sigma_eta_23;
lookup_tables.Sigma_eta_33 = Sigma_eta_33;
lookup_tables.mu_X_1 = mu_X_1;
lookup_tables.mu_X_2 = mu_X_2;
lookup_tables.Sigma_X_11 = Sigma_X_11;
lookup_tables.Sigma_X_12 = Sigma_X_12;
lookup_tables.Sigma_X_22 = Sigma_X_22;
lookup_tables.mu_CA50 = mu_CA50;
lookup_tables.Sigma_CA50 = Sigma_CA50;
lookup_tables.DI_QTY_interp = DI_QTY_interp;
lookup_tables.DI_SOI_interp = DI_SOI_interp;

end