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
    
    title(filename)
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

%% Joint Gaussian PDF

function plot_gaussian_pdf(Q_gross, X_res, mu, Sigma)

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

end

%% Start of Program

clear all
close all
clc

% myDir = '/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/';

% Constants
Pmax = 11; % MPa
gamma = 1.3;
EVC = -355;
EVO = 165;
CA_deg = linspace(-359.8, 360, 3600);
k = 1;
Q_min = 10000000; % Axes limits
Q_max = 0;
X_min = 10000000;
X_max = 0;

% Valve timings
i_evc = find(CA_deg == -355);
i_evo = find(CA_deg == 160.8000);

%% Loop through all files in the folder
for j = 1%1:length(myFiles)

    clear Pcyl_CA;
    clear n_cycles;

    % Load the file
    filename = [myDir, myFiles(j).name];
    load(filename)

    % extract the filename for output data
    file = erase(filename, ".mat");
    ofile = append(file, "_X_res.csv");
    fname = strrep(filename, '_', ' '); % replace _ with spaces for plot naming
    ofname = extractAfter(fname, '20'); % throw away everything before the date
    fname = append('20', ofname); % put the '20' back into '2024'

    n_cycles = length(Cylinder_1_Cycle_Data.IMEPn);

    % Plot the Cylinder Pressure vs. Crank Angle
   % figure(k); clf
   % plot_Pcyl(Pcyl_CA, Pmax, EVC, EVO);
   % k = k + 1;
   % 3600 data points per cycle = 7,200,000

    % Valve timings
    CoV_IMEP = std(Cylinder_1_Cycle_Data.IMEPn)/mean(Cylinder_1_Cycle_Data.IMEPn) * 100;
   
    % This fixes the reshape error for files with a different n_cycles
    if (mod(n_cycles, 100) == 0)
        % fprintf("%s\nCoV of IMEP: %.2f\n\n", filename, CoV_IMEP);
        % Only look at desirable conditions
        if (Cylinder_1_Cycle_Data.Injection_1_SOI(1) <= -40)
            fprintf("%s\n", filename);
            Pcyl_CA = reshape(Cylinder_1_Synch_Data.Cylinder_Pressure, [3600, n_cycles]) * 10 ^ -3;

            % Calculate the X_res (-)
            X_res = zeros([1, n_cycles]);
            for i = 1:n_cycles
                X_res(i) =  (Volume.Volume(i_evc) / Volume.Volume(i_evo)) * (Pcyl_CA(i_evc, i) / Pcyl_CA(i_evo, i)) ^ (1 / gamma);
            end
            
            % Write calculated values to a new file
            % writematrix(X_res, ofile);
            
            % Gross heat release (J)
            Q_gross = Cylinder_1_Cycle_Data.Gross_Heat_Release;

            % Residual gas fraction in percentage (%)
            X_res_per = X_res * 100;

            % Combine data into a single matrix
            data = [X_res_per; Q_gross]';

            % Estimate the mean and covariance matrix
            mu = mean(data);
            Sigma = cov(data);

            % Joint Gaussian PDF
            % X | Q = a ~ N(cond_mean, cond_variance)

            % Conditional variance
            cond_covariance = Sigma(1, 1) - Sigma(1, 2) * Sigma(2, 2)^-1 * Sigma(2, 1);

            % Known value of Q in the conditional distribution
            a = Q_gross;

            % Gaussian Conditional mean
            cond_mean = mu(1) + Sigma(1, 2) * Sigma(2, 2)^-1 * (a - mu(2));

            % Simulate residual gas fraction in percentage
            X_res_per_sim = normrnd(cond_mean, cond_covariance);

            % Kullback Liebler Divergence
            KL_div = 0;
            for i = 1:n_cycles
                KL_div = KL_div + X_res_per(i) * log(X_res_per(i) / X_res_per_sim(i));
            end

            fprintf("K-L Divergence Factor: %f\n\n", KL_div);

            % Plotting Comparison
            figure(k)
            scatter(Q_gross,X_res_per); hold on
            axis([190 1080 0 10])
            scatter(Q_gross, X_res_per_sim); legend('experiment', 'simulation')
            title(fname)
            xlabel('Q_{Gross} (J)'); ylabel('X_{res} (%)')
            figfile = append(file, "_scatter.jpg");
            saveas(figure(k), figfile);
            hold off
            k = k + 1;

            % Q-Q plot
            % figure(k); hold on
            % qqplot(X_res, X_res_per_sim);
            % title(fname)
            % qqfile = append(file, "_qqplot.jpg");
            % saveas(figure(k), qqfile)
            % k = k + 1;
            % hold off
            % 
            % % Histogram
            % figure(k); hold on
            % histogram(X_res_per);
            % histogram(X_res_per_sim); legend('experiment', 'simulation')
            % title(fname)
            % histofile = append(file, "_histogram.jpg");
            % saveas(figure(k), histofile)
            % hold off
            % k = k + 1;
           
        end
    else 
        fprintf("Bad file: %s\n# of Cycles: %f\n\n", filename, n_cycles);
    end
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