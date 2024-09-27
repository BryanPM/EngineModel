%% Tabula rasa
clear all
clc

%% Load data
myDir = 'Model_Data/';
myFiles_mat = dir(fullfile(myDir, 'DI*SOI*.mat'));
myFiles_csv = dir(fullfile(myDir, 'DI*SOI*.csv'));

% Initialize variables
mu_eta_1 = [];
mu_eta_2 = [];
mu_eta_3 = [];
Sigma_eta_11 = [];
Sigma_eta_12 = [];
Sigma_eta_13 = [];
Sigma_eta_22 = [];
Sigma_eta_23 = [];
Sigma_eta_33 = [];

mu_X_1 = [];
mu_X_2 = [];
Sigma_X_11 = [];
Sigma_X_12 = [];
Sigma_X_22 = [];

M_fuel_all = [];
M_air_all = [];
Q_gross_all = [];
m_diesel_all = [];
m_ammonia_all = [];
m_air_all = [];
CA50_all = [];
DI_SOI_all = [];
DI_DOI_all = [];

m_diesel_mean = [];
CA50_mean = [];
CA50_std = [];
DI_SOI_mean = [];
DI_DOI_mean = [];

DI_SOI_des = [];
DI_qty_des = [];

for i = 1:length(myFiles_mat)

% Load mat files
fileName_mat = myFiles_mat(i).name;
load([myDir, fileName_mat]);

% Extract parameters
mu_eta_1 = [mu_eta_1; eta_c_mu(1)];
mu_eta_2 = [mu_eta_2; eta_c_mu(2)];
mu_eta_3 = [mu_eta_3; eta_c_mu(3)];
Sigma_eta_11 = [Sigma_eta_11; eta_c_Sigma(1,1)];
Sigma_eta_12 = [Sigma_eta_12; eta_c_Sigma(1,2)];
Sigma_eta_13 = [Sigma_eta_13; eta_c_Sigma(1,3)];
Sigma_eta_22 = [Sigma_eta_22; eta_c_Sigma(2,2)];
Sigma_eta_23 = [Sigma_eta_23; eta_c_Sigma(2,3)];
Sigma_eta_33 = [Sigma_eta_33; eta_c_Sigma(3,3)];

mu_X_1 = [mu_X_1; X_res_mu(1)];
mu_X_2 = [mu_X_2; X_res_mu(2)];
Sigma_X_11 = [Sigma_X_11; X_res_Sigma(1,1)];
Sigma_X_12 = [Sigma_X_12; X_res_Sigma(1,2)];
Sigma_X_22 = [Sigma_X_22; X_res_Sigma(2,2)];

% Load csv files
fileName_csv = myFiles_csv(i).name;
cycleData = readmatrix([myDir, fileName_csv]);

% Extract time series
M_fuel_all = [M_fuel_all; cycleData(:,1)];
M_air_all = [M_air_all; cycleData(:,2)];
Q_gross_all = [Q_gross_all; cycleData(:,4)];
m_diesel_all = [m_diesel_all; cycleData(:,6)];
m_ammonia_all = [m_ammonia_all; cycleData(:,7)];
m_air_all = [m_air_all; cycleData(:,8)];
CA50_all = [CA50_all; cycleData(:,9)];
DI_SOI_all = [DI_SOI_all; cycleData(:,10)];
DI_DOI_all = [DI_DOI_all; cycleData(:,11)];

% Extract averages
m_diesel_mean = [m_diesel_mean; mean(cycleData(:,6))];
CA50_mean = [CA50_mean; mean(cycleData(:,9))];
CA50_std = [CA50_std; std(cycleData(:,9))];
DI_SOI_mean = [DI_SOI_mean; mean(cycleData(:,10))];
DI_DOI_mean = [DI_DOI_mean; mean(cycleData(:,11))];

% Extract inputs
DI_SOI_des = [DI_SOI_des; DI_timing_des];
DI_qty_des = [DI_qty_des; DI_quantity_des];

end

%% Regress model parameters

VARIABLE = {mu_eta_1; mu_eta_2; mu_eta_3; Sigma_eta_11; Sigma_eta_12;
            Sigma_eta_13; Sigma_eta_22; Sigma_eta_23; Sigma_eta_33;
            mu_X_1; mu_X_2; Sigma_X_11; Sigma_X_12; Sigma_X_22; CA50_mean;
            CA50_std};
VARNAME = {'mu_eta_1'; 'mu_eta_2'; 'mu_eta_3'; 'Sigma_eta_11'; 'Sigma_eta_12';
           'Sigma_eta_13'; 'Sigma_eta_22'; 'Sigma_eta_23'; 'Sigma_eta_33';
           'mu_X_1'; 'mu_X_2'; 'Sigma_X_11'; 'Sigma_X_12'; 'Sigma_X_22';
           'mu_CA50'; 'Sigma_CA50'};
VARTITLE = {'E[\eta_c] (%)'; 'E[M_{fuel}] (mg)'; 'E[M_{air}] (mg)';
            'Var[\eta_c]'; 'Cov[\eta_c, M_{fuel}]'; 'Cov[\eta_c, M_{air}]';
            'Var[M_{fuel}]'; 'Cov[M_{fuel}, M_{air}]'; 'Var[M_{air}]';
            'E[X_{res}] (%)'; 'E[Q_{gross}] (J)'; 'Var[X_{res}]';
            'Cov[X_{res}, Q_{gross}]'; 'Var[Q_{gross}]'; 'E[CA50] (aTDC)';
            'Var[CA50]'};

% Refine sampling
n_points = 100;
% [X, Y] = meshgrid(linspace(0.9*min(DI_qty_des), 1.1*max(DI_qty_des), n_points), ...
%                  linspace(0.9*min(DI_SOI_des), 1.1*max(DI_SOI_des), n_points+1));
[X, Y] = meshgrid(linspace(0.85*min(m_diesel_mean)*1e6, 1.15*max(m_diesel_mean)*1e6, n_points), ...
                  linspace(0.85*min(DI_SOI_mean), 1.15*max(DI_SOI_mean), n_points+1));
writematrix(X, 'Model_Data/DI_QTY_interp.csv');
writematrix(Y, 'Model_Data/DI_SOI_interp.csv');

for j = 1:numel(VARIABLE)

% Gather data
% [xData, yData, zData] = prepareSurfaceData(DI_qty_des, DI_SOI_des, VARIABLE{j});
[xData, yData, zData] = prepareSurfaceData(m_diesel_mean*1e6, DI_SOI_mean, VARIABLE{j});

% Set up fittype and options.
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.ExtrapolationMethod = 'none';
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

figure;
subplot( 2, 1, 1 );
plot( fitresult, [xData, yData], zData, 'Style', 'Contour' );
xlabel('Diesel Injection Quantity (mg)');
ylabel('Diesel Injection timing (deg bTDC)');
title(VARTITLE{j}); grid off; colorbar

subplot( 2, 1, 2 );
plot( fitresult, [xData, yData], zData );
xlabel('Injection Quantity (mg)');
ylabel('Injection timing (deg bTDC)');
zlabel(VARTITLE{j}); grid off

print(['Model_Plots/', VARNAME{j}], '-dpng', '-r300');

% Refine grid
Z = feval(fitresult, X, Y);

% Fit smoother model
ft = 'linearinterp';
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.ExtrapolationMethod = 'nearest';
opts.Normalize = 'on';

% Fit model to data.
[xData, yData, zData] = prepareSurfaceData( X, Y, Z );
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Save lookup matrix
Z = feval(fitresult, X, Y);
writematrix(Z, ['Model_Data/',  VARNAME{j}, '.csv'])

end

%% Meta-model

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
M_fuel_sim(1) = M_fuel_all(1);
M_air_sim(1) = M_air_all(1);

% Import lookup tables
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

for i = 1:n_cycles

    % Diesel inputs
    Diesel_SOI = DI_SOI_all(i);
    Diesel_fuel = m_diesel_all(i);

    % Other inputs
    Ammonia_fuel = m_ammonia_all(i);
    Fresh_air = m_air_all(i);
    
    % State
    state = [M_fuel_sim(i); M_air_sim(i)];
    
    % Joint mean for combustion efficiency distribution
    mu_eta_1_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_1, Diesel_fuel*1e6, Diesel_SOI);
    mu_eta_2_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_2, Diesel_fuel*1e6, Diesel_SOI);
    mu_eta_3_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_3, Diesel_fuel*1e6, Diesel_SOI);
    
    eta_c_mu = [mu_eta_1_eval, mu_eta_2_eval, mu_eta_3_eval];

    % Covariance matrix
    Sigma_eta_11_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_11, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_12_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_12, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_13_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_13, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_22_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_22, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_23_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_23, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_33_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_33, Diesel_fuel*1e6, Diesel_SOI);

    eta_c_Sigma = [Sigma_eta_11_eval, Sigma_eta_12_eval, Sigma_eta_13_eval;
                   Sigma_eta_12_eval, Sigma_eta_22_eval, Sigma_eta_23_eval;
                   Sigma_eta_13_eval, Sigma_eta_23_eval, Sigma_eta_33_eval];

    % Combustion efficiency
    eta_c_sim(i) = conditional_Gauss(eta_c_mu, eta_c_Sigma, state'*1e6) / 100;
    
    % Effective LHV
    Q_LHV_eff = (Diesel_fuel * Q_LHV_diesel + Ammonia_fuel * Q_LHV_ammonia) ./ (Diesel_fuel + Ammonia_fuel);

    % Gross heat release
    Q_gross_sim(i) = eta_c_sim(i) * M_fuel_sim(i) * Q_LHV_eff;

    % Mean for CA50 distribution
    mu_CA50_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_CA50, Diesel_fuel*1e6, Diesel_SOI);
    % Covariance matrix
    Sigma_CA50_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50, Diesel_fuel*1e6, Diesel_SOI);

    CA50_sim(i) = normrnd(mu_CA50_eval, Sigma_CA50_eval);

    % Joint mean for residual gas fraction
    mu_X_1_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_X_1, Diesel_fuel*1e6, Diesel_SOI);
    mu_X_2_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_X_2, Diesel_fuel*1e6, Diesel_SOI);

    X_res_mu = [mu_X_1_eval, mu_X_2_eval];

    % Covariance matrix
    Sigma_X_11_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_X_11, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_X_12_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_X_12, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_X_22_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_X_22, Diesel_fuel*1e6, Diesel_SOI);

    X_res_Sigma = [Sigma_X_11_eval, Sigma_X_12_eval;
                   Sigma_X_12_eval, Sigma_X_22_eval];

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
figure; clf; hold on
plot(Q_gross_all)
plot(Q_gross_sim)
ylabel('Q_{gross} (J)'); legend('Experimental', 'Simulated')
screenSize = get(0, 'ScreenSize');
originalHeight = 400; % You can adjust this height value as needed
set(gcf, 'Units', 'pixels', 'Position', [0, 0, screenSize(3), originalHeight]);

print('Model_Plots/Q_gross_comparison', '-dpng', '-r300');

figure; clf; hold on
plot(CA50_all)
plot(CA50_sim)
ylabel('CA50 (aTDC)'); legend('Experimental', 'Simulated')
screenSize = get(0, 'ScreenSize');
originalHeight = 400; % You can adjust this height value as needed
set(gcf, 'Units', 'pixels', 'Position', [0, 0, screenSize(3), originalHeight]);

print('Model_Plots/CA50_comparison', '-dpng', '-r300');

%%
function x1_sim = conditional_Gauss(mu, Sigma, x2)

% Conditional Gaussian
% x1 | x2 ~ N(mu_cond, sigma_cond)

% Conditional standard deviation
sigma_cond = sqrt(Sigma(1, 1) - Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * Sigma(2:end, 1));
if sigma_cond < 0
    sigma_cond = 0
end
% Gaussian Conditional mean
mu_cond = mu(1) + Sigma(1, 2:end) * Sigma(2:end, 2:end)^-1 * (x2 - mu(2:end))';

% Simulate residual gas fraction in percentage
x1_sim = normrnd(mu_cond, sigma_cond);
end