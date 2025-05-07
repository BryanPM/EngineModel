import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import norm

def load_experiments(myFiles, n_cycles_file):
    
    # Initialize measured variables
    M_diesel_exp = []
    M_NH3_exp = []
    M_air_exp = []
    m_diesel_exp = []
    m_NH3_exp = []
    m_air_exp = []
    Q_gross_exp = []
    CA50_exp = []
    DI_SOI_exp = []
    X_res_exp = []

    for i, fileName in enumerate(myFiles):
        cycleData = pd.read_csv(os.path.join(myDir, fileName)).values

        # Extract time series
        M_diesel_exp.extend(cycleData[:n_cycles_file, 0])
        M_NH3_exp.extend(cycleData[:n_cycles_file, 1])
        M_air_exp.extend(cycleData[:n_cycles_file, 2])
        m_diesel_exp.extend(cycleData[:n_cycles_file, 3])
        m_NH3_exp.extend(cycleData[:n_cycles_file, 4])
        m_air_exp.extend(cycleData[:n_cycles_file, 5])
        Q_gross_exp.extend(cycleData[:n_cycles_file, 6])
        CA50_exp.extend(cycleData[:n_cycles_file, 7])
        DI_SOI_exp.extend(cycleData[:n_cycles_file, 8])
        X_res_exp.extend(cycleData[:n_cycles_file, 9])

    # Convert to numpy arrays
    M_diesel_exp = np.array(M_diesel_exp)
    M_NH3_exp = np.array(M_NH3_exp)
    M_air_exp = np.array(M_air_exp)
    m_diesel_exp = np.array(m_diesel_exp)
    m_NH3_exp = np.array(m_NH3_exp)
    m_air_exp = np.array(m_air_exp)
    Q_gross_exp = np.array(Q_gross_exp)
    CA50_exp = np.array(CA50_exp)
    DI_SOI_exp = np.array(DI_SOI_exp)
    X_res_exp = np.array(X_res_exp)
    
    return M_diesel_exp, M_NH3_exp, M_air_exp, \
           m_diesel_exp, m_NH3_exp, m_air_exp, \
           Q_gross_exp, CA50_exp, DI_SOI_exp, X_res_exp
           
def load_lookup_tables(directory):
    # Load lookup tables from files
    lookup_tables = {}
    files = [f for f in os.listdir(directory) if (f.startswith('mu_') or f.startswith('Sigma_') or f.startswith('DI_'))]

    for file in files:
        key = file.replace('.csv', '')
        lookup_tables[key] = pd.read_csv(os.path.join(directory, file)).values

    return lookup_tables

def Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables):

    # Extract grid points
    DI_QTY_interp = lookup_tables['DI_QTY_interp'][0, :]
    DI_SOI_interp = lookup_tables['DI_SOI_interp'][:, 0]

     # Create the interpolators for each lookup table
    def interpolate(table):
        interp_func = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), table, method='linear', bounds_error=True, fill_value=None)
        return interp_func([[Diesel_fuel * 1e6, Diesel_SOI]])[0]  # Ensure correct input format

    eta_c_mu = np.array([interpolate(lookup_tables[f'mu_eta_{i}'].T) for i in range(1, 5)])
    eta_c_Sigma = np.array([[interpolate(lookup_tables[f'Sigma_eta_{i}{j}'].T) for j in range(1, 5)] for i in range(1, 5)])

    X_res_mu = np.array([interpolate(lookup_tables[f'mu_X_{i}'].T) for i in range(1, 3)])
    X_res_Sigma = np.array([[interpolate(lookup_tables[f'Sigma_X_{i}{j}'].T) for j in range(1, 3)] for i in range(1, 3)])

    CA50_mu = np.array([interpolate(lookup_tables[f'mu_CA50_{i}'].T) for i in range(1, 6)])
    CA50_Sigma = np.array([[interpolate(lookup_tables[f'Sigma_CA50_{i}{j}'].T) for j in range(1, 6)] for i in range(1, 6)])

    return eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, CA50_mu, CA50_Sigma

def conditional_Gauss(mu, Sigma, x2):
    # Conditional variance
    sigma_cond = np.sqrt(Sigma[0, 0] - Sigma[0, 1:] @ np.linalg.inv(Sigma[1:, 1:]) @ Sigma[1:, 0])
    sigma_cond = max(sigma_cond, 0)
    
    # Conditional mean
    mu_cond = mu[0] + Sigma[0, 1:] @ np.linalg.inv(Sigma[1:, 1:]) @ (x2 - mu[1:])
    
    return norm.rvs(mu_cond, sigma_cond)

# Data directory
myDir = 'Dual_Fuel_Data/'
myFiles = [f for f in os.listdir(myDir) if f.endswith('_Data.csv')]
myFiles.sort()

# Import lookup tables
lookup_tables = load_lookup_tables(myDir)

# Cycles per file
n_cycles_file = 100

# Load experimental data
M_diesel_exp, M_NH3_exp, M_air_exp, \
m_diesel_exp, m_NH3_exp, m_air_exp, \
Q_gross_exp, CA50_exp, DI_SOI_exp, X_res_exp = load_experiments(myFiles, n_cycles_file)

# Total number of cycles
n_cycles = len(Q_gross_exp)

# Initialize simulation variables
M_diesel_sim = np.zeros(n_cycles)
M_NH3_sim = np.zeros(n_cycles)
M_air_sim = np.zeros(n_cycles)
eta_c_sim = np.zeros(n_cycles)
X_res_sim = np.zeros(n_cycles)
Q_gross_sim = np.zeros(n_cycles)
CA50_sim = np.zeros(n_cycles)

# Initial condition
M_diesel_sim[0] = M_diesel_exp[0]
M_NH3_sim[0] = M_NH3_exp[0]
M_air_sim[0] = M_air_exp[0]

# Cost function
cost_c = np.zeros(n_cycles)

# Feedback from EONS
m_diesel_EONS = np.zeros(n_cycles)
SOI_diesel_EONS = np.zeros(n_cycles)
# State to EONS
M_diesel_EONS = np.zeros(n_cycles)
M_NH3_EONS = np.zeros(n_cycles)
M_air_EONS = np.zeros(n_cycles)
CA50_EONS = np.zeros(n_cycles)

# Model parameters
AFR_diesel = 14.5
AFR_NH3 = 6.04
Q_LHV_diesel = 44.1e6
Q_LHV_NH3 = 18.6e6

# Cost function parameters (NEED TO BE CALIBRATED)
alpha = np.mean(1 / M_diesel_exp)
beta = np.mean(1 / (1 - Q_gross_exp / (Q_LHV_diesel * M_diesel_exp + Q_LHV_NH3 * M_NH3_exp)))
gamma = 1 / (np.mean(7 - CA50_exp) ** 2)

for i in range(n_cycles):
    
    # ORCAS feedforward inputs
    m_diesel_FF = m_diesel_exp[i]
    m_NH3_FF = m_NH3_exp[i]
    m_air_FF = m_air_exp[i]
    SOI_diesel_FF = DI_SOI_exp[i]
    
    # EONS feedback
    m_diesel_FB = m_diesel_EONS[i]
    SOI_diesel_FB = SOI_diesel_EONS[i]
    
    # FF + FB diesel strategy
    m_diesel_total = m_diesel_FF + m_diesel_FB
    SOI_diesel_total = SOI_diesel_FF + SOI_diesel_FB

    # State
    state = np.array([M_diesel_sim[i], M_NH3_sim[i], M_air_sim[i]])
    state_CA50 = np.array([M_diesel_sim[i] * 1e6, M_NH3_sim[i] * 1e6, M_air_sim[i] * 1e4,
                           SOI_diesel_total])

    # Conditional distribution parameters
    eta_c_mu, eta_c_Sigma, \
    X_res_mu, X_res_Sigma, \
    CA50_mu, CA50_Sigma = Gauss_parameters(m_diesel_total, SOI_diesel_total, lookup_tables)

    # Simulate CA50
    CA50_sim[i] = conditional_Gauss(CA50_mu, CA50_Sigma, state_CA50)

    # Combustion efficiency
    eta_c_sim[i] = conditional_Gauss(eta_c_mu, eta_c_Sigma, state * 1e6) / 100

    # Gross heat release
    Q_gross_sim[i] = eta_c_sim[i] * (Q_LHV_diesel * M_diesel_sim[i] + Q_LHV_NH3 * M_NH3_sim[i])

    # Residual gas fraction
    X_res_sim[i] = conditional_Gauss(X_res_mu, X_res_Sigma, Q_gross_sim[i]) / 100

    # Cost function # FITNESS FUNCTION TO EVALUATE SNN PERFORMANCE
    cost_c[i] = alpha * M_diesel_sim[i] + \
                beta * (1 - eta_c_sim[i]) + \
                gamma * (7 - CA50_sim[i]) ** 2
                
    # Residual mass matrix
    Matrix_res = X_res_sim[i] * np.array([
        [      1 - eta_c_sim[i],                 0,           0],
        [            0,                 1 - eta_c_sim[i],     0],
        [-AFR_diesel * eta_c_sim[i], -AFR_NH3 * eta_c_sim[i], 1]
    ])

    # Fresh fuel and air
    input = np.array([m_diesel_total, m_NH3_FF, m_air_FF])
    
    # States passed from ORCAS to EONS
    M_diesel_EONS[i] = (M_diesel_sim[i] - m_diesel_FF) * 1e6
    M_NH3_EONS[i] = (M_NH3_sim[i] - m_NH3_FF) * 1e5
    M_air_EONS[i] = (M_air_sim[i] - m_air_FF) * 1e4
    CA50_EONS[i] = (CA50_sim[i] - 7) * 0.2

    if i < n_cycles - 1:
        # Calculate next cycle
        next_state = Matrix_res @ state + input
        M_diesel_sim[i + 1] = next_state[0]
        M_NH3_sim[i + 1] = next_state[1]
        M_air_sim[i + 1] = next_state[2]

#%% Plot results
def plot_timeseries(sim_data, exp_data, ylabel):
    plt.figure()
    plt.plot(exp_data, label='Experimental')
    plt.plot(sim_data, label='Simulated')
    plt.ylabel(ylabel)
    plt.title(f"RMSE = {np.sqrt(np.mean((exp_data - sim_data)**2)):.2f}")
    plt.legend()
    plt.xlabel('Cycles')
    plt.show() 

plot_timeseries(M_diesel_sim * 1e6, M_diesel_exp * 1e6, 'M_diesel (g)')
plot_timeseries(M_NH3_sim * 1e6, M_NH3_exp * 1e6, 'M_NH3 (g)')
plot_timeseries(M_air_sim * 1e6, M_air_exp * 1e6, 'M_air (g)')
plot_timeseries(Q_gross_sim, Q_gross_exp, 'Q_gross (J)')
plot_timeseries(CA50_sim, CA50_exp, 'CA50 (aTDC)')
plot_timeseries(X_res_sim * 100, X_res_exp * 100, 'X_res (%)')

plt.figure()
plt.plot(eta_c_sim * 100)
plt.ylabel('eta_c (%)')
plt.xlabel('Cycles')
plt.show() 

plt.figure()
plt.plot(cost_c)
plt.ylabel('Cost function')
plt.xlabel('Cycles')
plt.show() 

plt.figure()
plt.plot(m_diesel_EONS * 1e6, label='m_diesel_EONS (g)')
plt.plot(SOI_diesel_EONS, label='SOI_diesel_EONS (deg)')
plt.ylabel('Feedback from EONS')
plt.legend()
plt.xlabel('Cycles')
plt.show() 

plt.figure()
plt.plot(M_diesel_EONS, label='state_diesel_EONS')
plt.plot(M_NH3_EONS, label='state_NH3_EONS')
plt.plot(M_air_EONS, label='state_air_EONS')
plt.plot(CA50_EONS, label='state_air_EONS')
plt.ylabel('Feedback from EONS')
plt.legend()
plt.xlabel('Cycles')
plt.show() 