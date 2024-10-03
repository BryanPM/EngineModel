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


def load_experiments(myDir, n_cycles_file):

    # Initialize variables
    Q_gross_exp = []
    m_diesel_exp = []
    CA50_exp = []
    DI_SOI_exp = []

    # Path to measured data
    myFiles = [f for f in os.listdir(myDir) if f.endswith('_Data.csv')]

    # Load CSV files
    for fileName in myFiles:

        # Read CSV file
        cycleData = pd.read_csv(os.path.join(myDir, fileName)).values

        # Get initial condition (Kg)
        if len(Q_gross_exp) == 0:
            M_fuel_init = cycleData[0, 0]
            M_air_init = cycleData[1, 1]

        # Extract time series
        Q_gross_exp.append(cycleData[:n_cycles_file, 3])
        m_diesel_exp.append(cycleData[:n_cycles_file, 5])
        CA50_exp.append(cycleData[:n_cycles_file, 8])
        DI_SOI_exp.append(cycleData[:n_cycles_file, 9])

    # Convert to NumPy arrays
    Q_gross_exp = np.concatenate(Q_gross_exp)
    CA50_exp = np.concatenate(CA50_exp)
    DI_SOI_exp = np.concatenate(DI_SOI_exp)
    m_diesel_exp = np.concatenate(m_diesel_exp)
    
    return M_fuel_init, M_air_init, Q_gross_exp, CA50_exp, DI_SOI_exp, m_diesel_exp


def Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables):
    
    # Extract grid points
    DI_QTY_interp = lookup_tables['DI_QTY_interp'][0, :]
    DI_SOI_interp = lookup_tables['DI_SOI_interp'][:, 0]
    
    # Create the interpolators using RegularGridInterpolator for each lookup table
    interp_mu_eta_1 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_1'].T)
    interp_mu_eta_2 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_2'].T)
    interp_mu_eta_3 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_3'].T)
    
    interp_Sigma_eta_11 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_11'].T)
    interp_Sigma_eta_12 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_12'].T)
    interp_Sigma_eta_13 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_13'].T)
    interp_Sigma_eta_22 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_22'].T)
    interp_Sigma_eta_23 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_23'].T)
    interp_Sigma_eta_33 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_33'].T)
    
    interp_mu_X_1 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_X_1'].T)
    interp_mu_X_2 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_X_2'].T)
    
    interp_Sigma_X_11 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_X_11'].T)
    interp_Sigma_X_12 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_X_12'].T)
    interp_Sigma_X_22 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_X_22'].T)
    
    interp_mu_CA50 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50'].T)
    interp_Sigma_CA50 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50'].T)
    
    # Prepare point for interpolation
    point = [Diesel_fuel, Diesel_SOI]

    # Interpolate the values for the given Diesel_fuel and Diesel_SOI
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


# Path to lookup tables and experimental data
myDir = 'Model_Data/'

# Load lookup tables
lookup_tables = load_lookup_tables(myDir)

# Cycles per file
n_cycles_file = 100

# Load experimental data
M_fuel_init, M_air_init, Q_gross_exp, CA50_exp, DI_SOI_exp, m_diesel_exp = \
load_experiments(myDir, n_cycles_file)

# Simulator constants
Q_LHV_diesel = 44.1 * 1e6   # Lower heating value of Diesel (J)
Q_LHV_ammonia = 18.6 * 1e6  # Lower heating value of Ammonia (J)
n_cycles = len(Q_gross_exp) # Number of cycles in simulation
alpha = 0.02                # Cost function parameter

# Scale inputs to simulate how EONS will send messages to ORCAS
DI_SOI_EONS = (DI_SOI_exp - 42.56) / 6.72
m_diesel_EONS = (m_diesel_exp * 1e6 - 7.85) / 2.26

# Initialize variables
M_fuel_sim = np.zeros(n_cycles)
M_air_sim = np.zeros(n_cycles)
eta_c_sim = np.zeros(n_cycles)
X_res_sim = np.zeros(n_cycles)
Q_gross_sim = np.zeros(n_cycles)
CA50_sim = np.zeros(n_cycles)
cost_c = np.zeros(n_cycles)
cost_cummulative = 0
M_fuel_scaled = np.zeros(n_cycles)
M_air_scaled = np.zeros(n_cycles)
CA50_scaled = np.zeros(n_cycles)
Diesel_SOI_EONS = np.zeros(n_cycles)
Diesel_fuel_EONS = np.zeros(n_cycles)
Diesel_SOI = np.zeros(n_cycles)
Diesel_fuel = np.zeros(n_cycles)

# Initial condition
M_fuel_sim[0] = M_fuel_init
M_air_sim[0] = M_air_init

# Loop through cycles
for i in range(n_cycles):

    # Diesel feedback inputs: THESE ARE THE INPUTS YOU CAN MODIFY
    Diesel_SOI_EONS[i]  = DI_SOI_EONS[i]   # Keep between -1 and 1
    Diesel_fuel_EONS[i] = m_diesel_EONS[i] # Keep between -1 and 1

    # Diesel feedforward inputs: THESE WILL BE PROVIDED FROM ORCAS,USING CONSTANTS FOR NOW
    Diesel_SOI_FF = 42.56   # deg aTDC
    Diesel_fuel_FF = 7.85   # mg

    # Ammonia and air inputs: THESE WILL BE PROVIDED FROM ORCAS,USING CONSTANTS FOR NOW
    Ammonia_fuel_FF = 69    # mg
    Fresh_air_FF = 1600     # mg

    # Combine EONS feedback and feelforward diesel commands
    Diesel_SOI[i] = Diesel_SOI_EONS[i] * 6.72 + Diesel_SOI_FF
    Diesel_fuel[i] = Diesel_fuel_EONS[i] * 2.26 + Diesel_fuel_FF

    # Gaussian parameters
    eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, mu_CA50_eval, Sigma_CA50_eval \
        = Gauss_parameters(Diesel_fuel[i], Diesel_SOI[i], lookup_tables)
    
    # Mass state
    state = np.array([M_fuel_sim[i], M_air_sim[i]])
    # Combustion state
    CA50_sim[i] = np.random.normal(mu_CA50_eval, Sigma_CA50_eval)
    
    # Combustion efficiency
    eta_c_sim[i] = conditional_Gauss(eta_c_mu, eta_c_Sigma, state * 1e6) / 100

    # Effective LHV
    Q_LHV_eff = (Diesel_fuel[i] * Q_LHV_diesel + Ammonia_fuel_FF * Q_LHV_ammonia) / \
                (Diesel_fuel[i] + Ammonia_fuel_FF)

    # Gross heat release
    Q_gross_sim[i] = eta_c_sim[i] * M_fuel_sim[i] * Q_LHV_eff

    # Residual gas fraction
    X_res_sim[i] = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim[i]) / 100

    # Cost function # FITNESS FUNCTION TO EVALUATE SNN PERFORMANCE
    cost_c[i] = alpha * (CA50_sim[i] - 7.5) ** 2 + (eta_c_sim[i] - 1) ** 2
    cost_cummulative += cost_c[i]

    # Residual mass matrix
    Matrix_res = X_res_sim[i] * np.array([[1 - eta_c_sim[i], 0], [eta_c_sim[i], 1]])

    # Fresh fuel and air (Kg)
    input_mass = np.array([Diesel_fuel[i] + Ammonia_fuel_FF, Fresh_air_FF]) * 1e-6

    # Scale state: THIS ARE THE STATES YOU CAN USE FOR TRAINING THE SNN
    M_fuel_scaled[i] = M_fuel_sim[i] * 1e5 - 7.85  # Should be between -1 and 1
    M_air_scaled[i] = M_air_sim[i] *  1e4 - 16.6   # Should be between -1 and 1
    CA50_scaled[i] = CA50_sim[i] * 1e-1 - 0.61     # Should be between -1 and 1

    # Dynamics
    if i < n_cycles - 1:
        next_state = Matrix_res @ state + input_mass
        M_fuel_sim[i+1] = next_state[0]
        M_air_sim[i+1] = next_state[1]

# RMSE Calculation
RMSE_Q_gross = np.sqrt(np.mean((Q_gross_exp - Q_gross_sim) ** 2))
RMSE_CA50 = np.sqrt(np.mean((CA50_exp - CA50_sim) ** 2))

# Plot Q_gross
plt.figure()
plt.plot(Q_gross_exp, label='Experimental')
plt.plot(Q_gross_sim, label='Simulated')
plt.ylabel('Q_gross (J)')
plt.legend()
plt.title(f'Q_gross: RMSE = {RMSE_Q_gross:.1f} J')
plt.xlabel('Cycles')
plt.show()

# Plot CA50
plt.figure()
plt.plot(CA50_exp, label='Experimental')
plt.plot(CA50_sim, label='Simulated')
plt.ylabel('CA50 (aTDC)')
plt.legend()
plt.title(f'CA50: RMSE = {RMSE_CA50:.2f} deg')
plt.xlabel('Cycles')
plt.show()

# Plot combustion efficiency
plt.figure()
plt.plot(eta_c_sim)
plt.legend()
plt.title('Combustion efficiency')
plt.xlabel('Cycles')
plt.ylim([0, 1])
plt.show()

# Plot cost function
plt.figure()
plt.plot(cost_c)
plt.legend()
plt.title(f'Cost function: Total = {cost_cummulative:.1f}')
plt.xlabel('Cycles')
plt.show()

# Plot feedback from EONS
fig, ax = plt.subplots(2, 1)
ax[0].plot(Diesel_SOI_EONS)
ax[0].set_title('Feedback from EONS')
ax[0].set_ylabel('Diesel SOI')
ax[1].plot(Diesel_fuel_EONS)
ax[1].set_ylabel('Diesel mass')
ax[1].set_xlabel('Cycle')