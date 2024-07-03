clear all
close all
clc

myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/';
% myDir = '~/Applications/UTORII_DATA/';
myFiles = dir(fullfile(myDir, '*DI*SOI*.mat'));

% Constants
Pmax = 11; % MPa
gamma = 1.3;
EVC = -355;
EVO = 165;
CA_deg = linspace(-359.8, 360, 3600);

% Valve timings
i_evc = find(CA_deg == -355);
i_evo = find(CA_deg == 160.8000);

k = 1; % figure counter

%% Loop through all files in the folder
for j = 15%1:length(myFiles)

    clear Pcyl_CA;

    % Load the file
    filename = [myDir, myFiles(j).name];
    load(filename)

    % extract the filename for output data
    file = erase(filename, ".mat");
    ofile = append(file, "_X_res.csv");

    n_cycles = length(Cylinder_1_Synch_Data.Cylinder_Pressure) / 3600;
    Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure, [3600, n_cycles]) * 10 ^ -3;

    % Plot the Cylinder Pressure vs. Crank Angle
   % figure(k); clf
   % plot_Pcyl(Pcyl_CA, Pmax, EVC, EVO);
   % k = k + 1;
   % 2000 cycles, 3600 data points per cycle = 7,200,000

    % Valve timings
    CoV_IMEP = std(Cylinder_1_Cycle_Data.IMEPn)/mean(Cylinder_1_Cycle_Data.IMEPn) * 100;
   
    % Only look at low variable conditions
    if (CoV_IMEP <= 3)
        fprintf("Writing file %s to %s\nCoV of IMEP: %.2f\n\n", filename, ofile, CoV_IMEP);

        % Calculate the X_res as a function of cycle, k
        X_res = zeros([2000, 1]);
        for i = 1:n_cycles
            X_res(i) =  (Volume.Volume(i_evc) / Volume.Volume(i_evo)) * (Pcyl_CA(i_evc, i) / Pcyl_CA(i_evo, i)) ^ (1 / gamma);
        end

        % writematrix(X_res, ofile);
        % Plot after main loop
    end
end

plot_X_res(X_res, Cylinder_1_Cycle_Data.Gross_Heat_Release, filename, k)

%% Joint Gaussian PDF

Q_gross = Cylinder_1_Cycle_Data.Gross_Heat_Release';

% Combine data into a single matrix
data = [Q_gross, X_res];

% Estimate the mean and covariance matrix
mu = mean(data);
Sigma = cov(data);

% Create a grid of (x, y) values
x = linspace(min(Q_gross), max(Q_gross), 100);
y = linspace(min(X_res), max(X_res), 100);
[X, Y] = meshgrid(x, y);

% Compute the bivariate Gaussian PDF
Z = mvnpdf([X(:) Y(:)], mu, Sigma);
Z = reshape(Z, length(x), length(y));

% Plot the level curves (contour plot)
figure; hold on
scatter(Q_gross, X_res)
contour(X, Y, Z);
xlabel('Q_{gross}');
ylabel('X_{res}');
title('Level Curves of Fitted Bivariate Gaussian Distribution');
grid on;

%% Function to plot Gross_Heat_Release v. X_res

function plot_X_res(X_res, Gross_Heat_Release, filename, k)

    % Polynomial fit
    polymodel = polyfit(Gross_Heat_Release, X_res, 5);
    x = 770:(200/2000):1030; % X values
    f1 = polyval(polymodel, x);

    % Plotting
    figure(k); clf
    plot(Gross_Heat_Release, X_res, 'o');

    hold on

    plot(x, f1, '-');
    fname = strrep(filename, '_', ' ');
    title(fname)
    xlabel("Gross Heat Release [J]");
    ylabel("X res [%]");
    hold off
end

%% Function to plot, lifted from ANALYSIS.m
function plot_Pcyl(Pcyl_CA, Pmax, EVC, EVO)
    % Crank angle undisvidual cycle
    CA_deg = linspace(-359.8, 360, 3600);

    % Isolate compression and power stroke
    idx_comb  = find( CA_deg >= -360 & CA_deg <= 360 );
    CAD_comb  = CA_deg(idx_comb);
    
    % Combustion part
    Pcyl_comb = Pcyl_CA(idx_comb,:);
    
    % Plot
    colormap(summer);

    plot(CAD_comb,Pcyl_comb);
    xline(EVC, 'r', {'EVC'});
    xline(EVO, 'r', {'EVO'});
    xlabel('Crank Angle (deg)');
    ylabel('Cylinder Pressure (MPa)'); set(gca,'xtick',-400:60:400);
    xlim([-400 400]); ylim([0 Pmax]);
    title("Cylinder Pressure");
end