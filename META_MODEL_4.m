%% Tabula rasa
clear all
clc

%% Load data
myDir = 'Dual_Fuel_Data/';
myFiles_mat = dir(fullfile(myDir, 'DI*SOI*.mat'));
myFiles_csv = dir(fullfile(myDir, 'DI*SOI*.csv'));

% Initialize variables
mu_eta_1 = [];
mu_eta_2 = [];
mu_eta_3 = [];
mu_eta_4 = [];
Sigma_eta_11 = [];
Sigma_eta_12 = [];
Sigma_eta_13 = [];
Sigma_eta_14 = [];
Sigma_eta_22 = [];
Sigma_eta_23 = [];
Sigma_eta_24 = [];
Sigma_eta_33 = [];
Sigma_eta_34 = [];
Sigma_eta_44 = [];

mu_X_1 = [];
mu_X_2 = [];
Sigma_X_11 = [];
Sigma_X_12 = [];
Sigma_X_22 = [];

mu_CA50_1 = [];
mu_CA50_2 = [];
mu_CA50_3 = [];
mu_CA50_4 = [];
mu_CA50_5 = [];
Sigma_CA50_11 = [];
Sigma_CA50_12 = [];
Sigma_CA50_13 = [];
Sigma_CA50_14 = [];
Sigma_CA50_15 = [];
Sigma_CA50_22 = [];
Sigma_CA50_23 = [];
Sigma_CA50_24 = [];
Sigma_CA50_25 = [];
Sigma_CA50_33 = [];
Sigma_CA50_34 = [];
Sigma_CA50_35 = [];
Sigma_CA50_44 = [];
Sigma_CA50_45 = [];
Sigma_CA50_55 = [];


M_diesel_all = [];
M_NH3_all = [];
M_air_all = [];
Q_gross_all = [];
m_diesel_all = [];
m_NH3_all = [];
m_air_all = [];
CA50_all = [];
DI_SOI_all = [];

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
mu_eta_4 = [mu_eta_4; eta_c_mu(4)];
Sigma_eta_11 = [Sigma_eta_11; eta_c_Sigma(1,1)];
Sigma_eta_12 = [Sigma_eta_12; eta_c_Sigma(1,2)];
Sigma_eta_13 = [Sigma_eta_13; eta_c_Sigma(1,3)];
Sigma_eta_14 = [Sigma_eta_14; eta_c_Sigma(1,4)];
Sigma_eta_22 = [Sigma_eta_22; eta_c_Sigma(2,2)];
Sigma_eta_23 = [Sigma_eta_23; eta_c_Sigma(2,3)];
Sigma_eta_24 = [Sigma_eta_24; eta_c_Sigma(2,4)];
Sigma_eta_33 = [Sigma_eta_33; eta_c_Sigma(3,3)];
Sigma_eta_34 = [Sigma_eta_34; eta_c_Sigma(3,4)];
Sigma_eta_44 = [Sigma_eta_44; eta_c_Sigma(4,4)];

mu_X_1 = [mu_X_1; X_res_mu(1)];
mu_X_2 = [mu_X_2; X_res_mu(2)];
Sigma_X_11 = [Sigma_X_11; X_res_Sigma(1,1)];
Sigma_X_12 = [Sigma_X_12; X_res_Sigma(1,2)];
Sigma_X_22 = [Sigma_X_22; X_res_Sigma(2,2)];

mu_CA50_1 = [mu_CA50_1; CA50_mu(1)];
mu_CA50_2 = [mu_CA50_2; CA50_mu(2)];
mu_CA50_3 = [mu_CA50_3; CA50_mu(3)];
mu_CA50_4 = [mu_CA50_4; CA50_mu(4)];
mu_CA50_5 = [mu_CA50_5; CA50_mu(5)];
Sigma_CA50_11 = [Sigma_CA50_11; CA50_Sigma(1,1)];
Sigma_CA50_12 = [Sigma_CA50_12; CA50_Sigma(1,2)];
Sigma_CA50_13 = [Sigma_CA50_13; CA50_Sigma(1,3)];
Sigma_CA50_14 = [Sigma_CA50_14; CA50_Sigma(1,4)];
Sigma_CA50_15 = [Sigma_CA50_15; CA50_Sigma(1,5)];
Sigma_CA50_22 = [Sigma_CA50_22; CA50_Sigma(2,2)];
Sigma_CA50_23 = [Sigma_CA50_23; CA50_Sigma(2,3)];
Sigma_CA50_24 = [Sigma_CA50_24; CA50_Sigma(2,4)];
Sigma_CA50_25 = [Sigma_CA50_25; CA50_Sigma(2,5)];
Sigma_CA50_33 = [Sigma_CA50_33; CA50_Sigma(3,3)];
Sigma_CA50_34 = [Sigma_CA50_34; CA50_Sigma(3,4)];
Sigma_CA50_35 = [Sigma_CA50_35; CA50_Sigma(3,5)];
Sigma_CA50_44 = [Sigma_CA50_44; CA50_Sigma(4,4)];
Sigma_CA50_45 = [Sigma_CA50_45; CA50_Sigma(4,5)];
Sigma_CA50_55 = [Sigma_CA50_55; CA50_Sigma(5,5)];

% Load csv files
fileName_csv = myFiles_csv(i).name;
cycleData = readmatrix([myDir, fileName_csv]);

% Extract time series
M_diesel_all = [M_diesel_all; cycleData(:,1)];
M_NH3_all = [M_NH3_all; cycleData(:,2)];
M_air_all = [M_air_all; cycleData(:,3)];
m_diesel_all = [m_diesel_all; cycleData(:,4)];
m_NH3_all = [m_NH3_all; cycleData(:,5)];
m_air_all = [m_air_all; cycleData(:,6)];

Q_gross_all = [Q_gross_all; cycleData(:,7)];
CA50_all = [CA50_all; cycleData(:,8)];
DI_SOI_all = [DI_SOI_all; cycleData(:,9)];

% Extract inputs
DI_SOI_des = [DI_SOI_des; DI_timing_des];
DI_qty_des = [DI_qty_des; DI_quantity_des];

end

%% Regress model parameters

VARIABLE = {mu_eta_1; mu_eta_2; mu_eta_3; mu_eta_4;
            Sigma_eta_11; Sigma_eta_12; Sigma_eta_13; Sigma_eta_14;
            Sigma_eta_22; Sigma_eta_23; Sigma_eta_24;
            Sigma_eta_33; Sigma_eta_34;
            Sigma_eta_44;
            mu_X_1; mu_X_2;
            Sigma_X_11; Sigma_X_12;
            Sigma_X_22;
            mu_CA50_1; mu_CA50_2; mu_CA50_3; mu_CA50_4; mu_CA50_5;
            Sigma_CA50_11; Sigma_CA50_12; Sigma_CA50_13; Sigma_CA50_14; Sigma_CA50_15;
            Sigma_CA50_22; Sigma_CA50_23; Sigma_CA50_24; Sigma_CA50_25;
            Sigma_CA50_33; Sigma_CA50_34; Sigma_CA50_35; 
            Sigma_CA50_44; Sigma_CA50_45;
            Sigma_CA50_55};
VARNAME = {'mu_eta_1'; 'mu_eta_2'; 'mu_eta_3'; 'mu_eta_4';
           'Sigma_eta_11'; 'Sigma_eta_12'; 'Sigma_eta_13'; 'Sigma_eta_14';
           'Sigma_eta_22'; 'Sigma_eta_23'; 'Sigma_eta_24';
           'Sigma_eta_33'; 'Sigma_eta_34';
           'Sigma_eta_44';
           'mu_X_1'; 'mu_X_2';
           'Sigma_X_11'; 'Sigma_X_12';
           'Sigma_X_22';
           'mu_CA50_1'; 'mu_CA50_2'; 'mu_CA50_3'; 'mu_CA50_4'; 'mu_CA50_5';
           'Sigma_CA50_11'; 'Sigma_CA50_12'; 'Sigma_CA50_13'; 'Sigma_CA50_14'; 'Sigma_CA50_15';
           'Sigma_CA50_22'; 'Sigma_CA50_23'; 'Sigma_CA50_24'; 'Sigma_CA50_25';
           'Sigma_CA50_33'; 'Sigma_CA50_34'; 'Sigma_CA50_35';
           'Sigma_CA50_44'; 'Sigma_CA50_45';
           'Sigma_CA50_55'};
VARTITLE = {'E[\eta_c] (%)'; 'E[M_{diesel}] (mg)'; 'E[M_{NH3}] (mg)'; 'E[M_{air}] (mg)';
            'Var[\eta_c]'; 'Cov[\eta_c, M_{diesel}]'; 'Cov[\eta_c, M_{NH3}]'; 'Cov[\eta_c, M_{air}]';
            'Var[M_{diesel}]'; 'Cov[M_{diesel}, M_{NH3}]'; 'Cov[M_{diesel}, M_{air}]';
            'Var[M_{NH3}]'; 'Cov[M_{NH3}, M_{air}]';
            'Var[M_{air}]';
            'E[X_{res}] (%)'; 'E[Q_{gross}] (J)'; 'Var[X_{res}]';
            'Cov[X_{res}, Q_{gross}]'; 'Var[Q_{gross}]';
            'E[CA50] (aTDC)'; 'E[M_{diesel}] (mg)'; 'E[M_{NH3}] (mg)'; 'E[M_{air}] (dg)'; 'E[SOI] (deg bTDC)'
            'Var[CA50]'; 'Cov[CA50, M_{diesel}]'; 'Cov[CA50, M_{NH3}]'; 'Cov[CA50, M_{air}]'; 'Cov[CA50, SOI]';
            'Var[M_{diesel}]'; 'Cov[M_{diesel}, M_{NH3}]'; 'Cov[M_{diesel}, M_{air}]'; 'Cov[M_{diesel}, SOI]';
            'Var[M_{NH3}]'; 'Cov[M_{NH3}, M_{air}]'; 'Cov[M_{NH3}, SOI]';
            'Var[M_{air}]'; 'Cov[M_{air}, SOI]';
            'Var[SOI]'};

% Refine sampling
n_points = 100;
[X, Y] = meshgrid(linspace(0.9*min(DI_qty_des), 1.1*max(DI_qty_des), n_points), ...
                  linspace(0.9*min(DI_SOI_des), 1.1*max(DI_SOI_des), n_points+1));
writematrix(X, 'Dual_Fuel_Data/DI_QTY_interp.csv');
writematrix(Y, 'Dual_Fuel_Data/DI_SOI_interp.csv');

for j = 1:numel(VARIABLE)

% Gather data
[xData, yData, zData] = prepareSurfaceData(DI_qty_des, DI_SOI_des, VARIABLE{j});

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

print(['Dual_Fuel_Plots/', VARNAME{j}], '-dpng', '-r300');

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
writematrix(Z, ['Dual_Fuel_Data/',  VARNAME{j}, '.csv'])

end

%% Meta-model

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
M_diesel_sim(1) = M_diesel_all(1);
M_NH3_sim(1) = M_NH3_all(1);
M_air_sim(1) = M_air_all(1);

% Import lookup tables
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

% Model parameters
AFR_diesel = 14.5;
AFR_NH3    = 6.04;
Q_LHV_diesel = 44.1e6;
Q_LHV_NH3    = 18.6e6;

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
    
    % Joint mean for combustion efficiency distribution
    mu_eta_1_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_1, Diesel_fuel*1e6, Diesel_SOI);
    mu_eta_2_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_2, Diesel_fuel*1e6, Diesel_SOI);
    mu_eta_3_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_3, Diesel_fuel*1e6, Diesel_SOI);
    mu_eta_4_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_eta_4, Diesel_fuel*1e6, Diesel_SOI);
    
    eta_c_mu = [mu_eta_1_eval, mu_eta_2_eval, mu_eta_3_eval, mu_eta_4_eval];

    % Covariance matrix
    Sigma_eta_11_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_11, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_12_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_12, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_13_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_13, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_14_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_14, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_22_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_22, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_23_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_23, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_24_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_24, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_33_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_33, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_34_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_34, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_eta_44_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_eta_44, Diesel_fuel*1e6, Diesel_SOI);

    eta_c_Sigma = [Sigma_eta_11_eval, Sigma_eta_12_eval, Sigma_eta_13_eval, Sigma_eta_14_eval;
                   Sigma_eta_12_eval, Sigma_eta_22_eval, Sigma_eta_23_eval, Sigma_eta_24_eval;
                   Sigma_eta_13_eval, Sigma_eta_23_eval, Sigma_eta_33_eval, Sigma_eta_34_eval;
                   Sigma_eta_14_eval, Sigma_eta_24_eval, Sigma_eta_34_eval, Sigma_eta_44_eval];

    % Joint mean for CA50 distribution
    mu_CA50_1_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_CA50_1, Diesel_fuel*1e6, Diesel_SOI);
    mu_CA50_2_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_CA50_2, Diesel_fuel*1e6, Diesel_SOI);
    mu_CA50_3_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_CA50_3, Diesel_fuel*1e6, Diesel_SOI);
    mu_CA50_4_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_CA50_4, Diesel_fuel*1e6, Diesel_SOI);
    mu_CA50_5_eval = interp2(DI_QTY_interp, DI_SOI_interp, mu_CA50_5, Diesel_fuel*1e6, Diesel_SOI);

    CA50_mu = [mu_CA50_1_eval, mu_CA50_2_eval, mu_CA50_3_eval, mu_CA50_4_eval, mu_CA50_5_eval];

    % Covariance matrix
    Sigma_CA50_11_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_11, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_12_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_12, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_13_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_13, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_14_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_14, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_15_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_15, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_22_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_22, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_23_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_23, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_24_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_24, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_25_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_25, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_33_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_33, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_34_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_34, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_35_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_35, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_44_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_44, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_45_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_45, Diesel_fuel*1e6, Diesel_SOI);
    Sigma_CA50_55_eval = interp2(DI_QTY_interp, DI_SOI_interp, Sigma_CA50_55, Diesel_fuel*1e6, Diesel_SOI);

    CA50_Sigma = [Sigma_CA50_11_eval, Sigma_CA50_12_eval, Sigma_CA50_13_eval, Sigma_CA50_14_eval, Sigma_CA50_15_eval;
                  Sigma_CA50_12_eval, Sigma_CA50_22_eval, Sigma_CA50_23_eval, Sigma_CA50_24_eval, Sigma_CA50_25_eval;
                  Sigma_CA50_13_eval, Sigma_CA50_23_eval, Sigma_CA50_33_eval, Sigma_CA50_34_eval, Sigma_CA50_35_eval;
                  Sigma_CA50_14_eval, Sigma_CA50_24_eval, Sigma_CA50_34_eval, Sigma_CA50_44_eval, Sigma_CA50_45_eval;
                  Sigma_CA50_15_eval, Sigma_CA50_25_eval, Sigma_CA50_35_eval, Sigma_CA50_45_eval, Sigma_CA50_55_eval];

    % CA50
    CA50_sim(i) = conditional_Gauss(CA50_mu, CA50_Sigma, state_CA50);

    % Combustion efficiency
    eta_c_sim(i) = conditional_Gauss(eta_c_mu, eta_c_Sigma, state'*1e6) / 100;
    
    % Gross heat release
    Q_gross_sim(i) = eta_c_sim(i) * (Q_LHV_diesel * M_diesel_sim(i) + Q_LHV_NH3 * M_NH3_sim(i));

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

%% Plot
RMSE = sqrt(mean((Q_gross_all - Q_gross_sim).^2));
figure; clf; hold on; box on;
plot(Q_gross_all)
plot(Q_gross_sim)
ylabel('Q_{gross} (J)'); legend('Experimental', 'Simulated')
title("RMSE = " + num2str(RMSE, 3) + " (J)")
screenSize = get(0, 'ScreenSize'); xlabel('Cycles');
originalHeight = 400; % You can adjust this height value as needed
set(gcf, 'Units', 'pixels', 'Position', [0, 0, screenSize(3), originalHeight]);

print('Dual_Fuel_Plots/Q_gross_comparison', '-dpng', '-r300');

RMSE = sqrt(mean((CA50_all - CA50_sim).^2));
figure; clf; hold on; box on;
plot(CA50_all)
plot(CA50_sim)
ylabel('CA50 (aTDC)'); legend('Experimental', 'Simulated')
title("RMSE = " + num2str(RMSE, 2) + " (deg)")
screenSize = get(0, 'ScreenSize'); xlabel('Cycles');
originalHeight = 400; % You can adjust this height value as needed
set(gcf, 'Units', 'pixels', 'Position', [0, 0, screenSize(3), originalHeight]);

print('Dual_Fuel_Plots/CA50_comparison', '-dpng', '-r300');

%%
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