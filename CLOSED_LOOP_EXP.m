load('/Users/1pq/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Research/NTRC/UTORII_2024/UTORII Data/20241028/20241028_NH3_DIQty8_SOI45_Network5_Switch_007.mat')

signal = Cylinder_1_Cycle_Data.IMEPg;

% Define moving window size
window_size = 200; % Window size for moving average

% Compute moving mean and moving standard deviation
moving_mean = movmean(signal, [window_size, 0], "Endpoints", "fill");
moving_std = movstd(signal, [window_size, 0]);

% Compute coefficient of variation (CV)
cv_values = moving_std ./ moving_mean * 100;

idx = find(sgolayfilt(Cylinder_1_Cycle_Data.Injection_1_Duration, 3, 21) > 0.42, 1);
%%
figure(1); clf

subplot(6,1,1)
plot(-Cylinder_1_Cycle_Data.Injection_1_SOI)
ylabel("SOI^{dsl} (deg)")
xlim([0 3000])
xticklabels([]);


subplot(6,1,2)
plot(Cylinder_1_Cycle_Data.Feedforward_Fuel + [zeros(1,idx), Cylinder_1_Cycle_Data.EONS_Fuel_Qty(idx+1:end)])
ylabel("m_{in}^{dsl} (mg)")
xlim([0 3000])
xticklabels([]);


subplot(6,1,3)
plot(Cylinder_1_Cycle_Data.CA50)
ylabel("CA50 (deg)")
xlim([0 3000])
xticklabels([]);


subplot(6,1,4)
plot(Cylinder_1_Cycle_Data.Combustion_Efficiency)
ylabel("\eta_c (%)")
xlim([0 3000])
xticklabels([]);


subplot(6,1,5)
plot(Cylinder_1_Cycle_Data.Cost_C)
ylabel("Cost")
xlim([0 3000])
xticklabels([]);


subplot(6,1,6)
plot(cv_values); hold on;
plot([1 3000], [3 3], '--')
legend('Moving average estimate', 'Limit')
ylabel("CoV IMEP (%)")
xlim([0 3000])
xlabel('Combustion Cycles')