import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

def load_lookup_tables(myDir):
    
    # Load lookup tables from files
    lookup_tables = {}
    lookup_tables['mu_eta_1'] = pd.read_csv(os.path.join(myDir, 'mu_eta_1.csv'), header=None).values
    lookup_tables['mu_eta_2'] = pd.read_csv(os.path.join(myDir, 'mu_eta_2.csv'), header=None).values
    lookup_tables['mu_eta_3'] = pd.read_csv(os.path.join(myDir, 'mu_eta_3.csv'), header=None).values
    lookup_tables['Sigma_eta_11'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_11.csv'), header=None).values
    lookup_tables['Sigma_eta_12'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_12.csv'), header=None).values
    lookup_tables['Sigma_eta_13'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_13.csv'), header=None).values
    lookup_tables['Sigma_eta_22'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_22.csv'), header=None).values
    lookup_tables['Sigma_eta_23'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_23.csv'), header=None).values
    lookup_tables['Sigma_eta_33'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_33.csv'), header=None).values
    lookup_tables['mu_X_1'] = pd.read_csv(os.path.join(myDir, 'mu_X_1.csv'), header=None).values
    lookup_tables['mu_X_2'] = pd.read_csv(os.path.join(myDir, 'mu_X_2.csv'), header=None).values
    lookup_tables['Sigma_X_11'] = pd.read_csv(os.path.join(myDir, 'Sigma_X_11.csv'), header=None).values
    lookup_tables['Sigma_X_12'] = pd.read_csv(os.path.join(myDir, 'Sigma_X_12.csv'), header=None).values
    lookup_tables['Sigma_X_22'] = pd.read_csv(os.path.join(myDir, 'Sigma_X_22.csv'), header=None).values
    lookup_tables['mu_CA50'] = pd.read_csv(os.path.join(myDir, 'mu_CA50.csv'), header=None).values
    lookup_tables['Sigma_CA50'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50.csv'), header=None).values
    lookup_tables['DI_QTY_interp'] = pd.read_csv(os.path.join(myDir, 'DI_QTY_interp.csv'), header=None).values
    lookup_tables['DI_SOI_interp'] = pd.read_csv(os.path.join(myDir, 'DI_SOI_interp.csv'), header=None).values

    return lookup_tables

import numpy as np
from scipy.interpolate import RegularGridInterpolator

def Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables):
    # Convert Diesel_fuel from kg to mg for interpolation
    Diesel_fuel_mg = Diesel_fuel * 1e6
    
    # Ensure that the grid points are sorted in ascending order
    DI_QTY_interp = lookup_tables['DI_QTY_interp']
    DI_SOI_interp = lookup_tables['DI_SOI_interp']
    
    # If DI_QTY_interp and DI_SOI_interp are already in sorted order,
    # you can skip sorting. Otherwise, you can sort them as needed.
    
    # Create the interpolators using RegularGridInterpolator for each lookup table
    interp_mu_eta_1 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['mu_eta_1'].T)
    interp_mu_eta_2 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['mu_eta_2'].T)
    interp_mu_eta_3 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['mu_eta_3'].T)
    
    # Create interpolators for covariance matrices similarly
    interp_Sigma_eta_11 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_eta_11'].T)
    interp_Sigma_eta_12 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_eta_12'].T)
    interp_Sigma_eta_13 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_eta_13'].T)
    interp_Sigma_eta_22 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_eta_22'].T)
    interp_Sigma_eta_23 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_eta_23'].T)
    interp_Sigma_eta_33 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_eta_33'].T)
    
    interp_mu_X_1 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['mu_X_1'].T)
    interp_mu_X_2 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['mu_X_2'].T)
    
    interp_Sigma_X_11 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_X_11'].T)
    interp_Sigma_X_12 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_X_12'].T)
    interp_Sigma_X_22 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_X_22'].T)
    
    interp_mu_CA50 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['mu_CA50'].T)
    interp_Sigma_CA50 = RegularGridInterpolator((DI_QTY_interp[0, :], DI_SOI_interp[:, 0]), lookup_tables['Sigma_CA50'].T)
    
    # Prepare point for interpolation
    point = np.ascontiguousarray([Diesel_fuel_mg, Diesel_SOI])  # Ensure point is C-contiguous

    # Interpolate the values for the given Diesel_fuel_mg and Diesel_SOI
    mu_eta_1_eval = interp_mu_eta_1(point)
    mu_eta_2_eval = interp_mu_eta_2(point)
    mu_eta_3_eval = interp_mu_eta_3(point)

    Sigma_eta_11_eval = interp_Sigma_eta_11(point)
    Sigma_eta_12_eval = interp_Sigma_eta_12(point)
    Sigma_eta_13_eval = interp_Sigma_eta_13(point)
    Sigma_eta_22_eval = interp_Sigma_eta_22(point)
    Sigma_eta_23_eval = interp_Sigma_eta_23(point)
    Sigma_eta_33_eval = interp_Sigma_eta_33(point)

    mu_X_1_eval = interp_mu_X_1(point)
    mu_X_2_eval = interp_mu_X_2(point)

    Sigma_X_11_eval = interp_Sigma_X_11(point)
    Sigma_X_12_eval = interp_Sigma_X_12(point)
    Sigma_X_22_eval = interp_Sigma_X_22(point)

    mu_CA50_eval = interp_mu_CA50(point)
    Sigma_CA50_eval = interp_Sigma_CA50(point)

    # Pack outputs
    X_res_mu = np.array([mu_X_1_eval, mu_X_2_eval]).flatten()
    X_res_Sigma = np.array([[Sigma_X_11_eval, Sigma_X_12_eval], [Sigma_X_12_eval, Sigma_X_22_eval]]).reshape(2, 2)

    eta_c_mu = np.array([mu_eta_1_eval, mu_eta_2_eval, mu_eta_3_eval]).flatten()
    eta_c_Sigma = np.array([[Sigma_eta_11_eval, Sigma_eta_12_eval, Sigma_eta_13_eval],
                            [Sigma_eta_12_eval, Sigma_eta_22_eval, Sigma_eta_23_eval],
                            [Sigma_eta_13_eval, Sigma_eta_23_eval, Sigma_eta_33_eval]]).reshape(3, 3)

    return eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, mu_CA50_eval.flatten()[0], Sigma_CA50_eval.flatten()[0]



def conditional_Gauss(mu, Sigma, x2):
    
    # Conditional variance
    sigma_cond = np.sqrt(Sigma[0, 0] - Sigma[0, 1:] @ np.linalg.inv(Sigma[1:, 1:]) @ Sigma[1:, 0])
    # Avoid negative variances
    sigma_cond = max(sigma_cond, 0)

    # Conditional mean
    mu_cond = mu[0] + Sigma[0, 1:] @ np.linalg.inv(Sigma[1:, 1:]) @ (x2 - mu[1:])

    return np.random.normal(mu_cond, sigma_cond)


# Path to measured data
myDir = 'Model_Data/'
myFiles = [f for f in os.listdir(myDir) if f.endswith('_Data.csv')]

# Initialize variables
Q_gross_all = []
m_diesel_all = []
m_ammonia_all = []
m_air_all = []
CA50_all = []
DI_SOI_all = []

# Cycles per file
n_cycles_file = 100

# Load CSV files
for fileName in myFiles:

    cycleData = pd.read_csv(os.path.join(myDir, fileName)).values

    # Initial condition
    if len(Q_gross_all) == 0:
        M_fuel_init = cycleData[0, 0]
        M_air_init = cycleData[1, 1]

    # Extract time series
    Q_gross_all.append(cycleData[:n_cycles_file, 3])
    m_diesel_all.append(cycleData[:n_cycles_file, 5])
    m_ammonia_all.append(cycleData[:n_cycles_file, 6])
    m_air_all.append(cycleData[:n_cycles_file, 7])
    CA50_all.append(cycleData[:n_cycles_file, 8])
    DI_SOI_all.append(cycleData[:n_cycles_file, 9])

# Convert to NumPy arrays
Q_gross_all = np.concatenate(Q_gross_all)
m_diesel_all = np.concatenate(m_diesel_all)
m_ammonia_all = np.concatenate(m_ammonia_all)
m_air_all = np.concatenate(m_air_all)
CA50_all = np.concatenate(CA50_all)
DI_SOI_all = np.concatenate(DI_SOI_all)

# Simulator constants
Q_LHV_diesel = 44.1 * 1e6
Q_LHV_ammonia = 18.6 * 1e6
n_cycles = len(Q_gross_all)

# Initialize variables
M_fuel_sim = np.zeros(n_cycles)
M_air_sim = np.zeros(n_cycles)
eta_c_sim = np.zeros(n_cycles)
X_res_sim = np.zeros(n_cycles)
Q_gross_sim = np.zeros(n_cycles)
CA50_sim = np.zeros(n_cycles)

# Initial condition
M_fuel_sim[0] = M_fuel_init
M_air_sim[0] = M_air_init

# Load lookup tables
lookup_tables = load_lookup_tables(myDir)

for i in range(n_cycles):

    # Diesel inputs: THIS ARE THE INPUTS YOU CAN MODIFY
    Diesel_SOI = DI_SOI_all[i]
    Diesel_fuel = m_diesel_all[i]

    # Ammonia and air inputs
    Ammonia_fuel = m_ammonia_all[i]
    Fresh_air = m_air_all[i]

    # State
    state = np.array([M_fuel_sim[i], M_air_sim[i]])

    # Gaussian parameters
    eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, mu_CA50_eval, Sigma_CA50_eval = Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables)

    # Simulate CA50
    CA50_sim[i] = np.random.normal(mu_CA50_eval, Sigma_CA50_eval)

    # Combustion efficiency
    eta_c_sim[i] = conditional_Gauss(eta_c_mu, eta_c_Sigma, state * 1e6) / 100

    # Effective LHV
    Q_LHV_eff = (Diesel_fuel * Q_LHV_diesel + Ammonia_fuel * Q_LHV_ammonia) / (Diesel_fuel + Ammonia_fuel)

    # Gross heat release
    Q_gross_sim[i] = eta_c_sim[i] * M_fuel_sim[i] * Q_LHV_eff

    # Residual gas fraction
    X_res_sim[i] = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim[i]) / 100

    # Residual mass matrix
    Matrix_res = X_res_sim[i] * np.array([[1 - eta_c_sim[i], 0], [eta_c_sim[i], 1]])

    # Fresh fuel and air
    input_state = np.array([Diesel_fuel + Ammonia_fuel, Fresh_air])

    if i < n_cycles - 1:
        next_state = Matrix_res @ state + input_state
        M_fuel_sim[i+1] = next_state[0]
        M_air_sim[i+1] = next_state[1]

# RMSE Calculation
RMSE_Q_gross = np.sqrt(np.mean((Q_gross_all - Q_gross_sim) ** 2))
RMSE_CA50 = np.sqrt(np.mean((CA50_all - CA50_sim) ** 2))

# Plot Q_gross
plt.figure()
plt.plot(Q_gross_all, label='Experimental')
plt.plot(Q_gross_sim, label='Simulated')
plt.ylabel('Q_gross (J)')
plt.legend()
plt.title(f'Q_gross: RMSE = {RMSE_Q_gross:.1f} J')
plt.xlabel('Cycles')
plt.show()

# Plot CA50
plt.figure()
plt.plot(CA50_all, label='Experimental')
plt.plot(CA50_sim, label='Simulated')
plt.ylabel('CA50 (aTDC)')
plt.legend()
plt.title(f'CA50: RMSE = {RMSE_CA50:.2f} deg')
plt.xlabel('Cycles')
plt.show()