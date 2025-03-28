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
    lookup_tables['mu_eta_4'] = pd.read_csv(os.path.join(myDir, 'mu_eta_4.csv'), header=None).values
    lookup_tables['Sigma_eta_11'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_11.csv'), header=None).values
    lookup_tables['Sigma_eta_12'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_12.csv'), header=None).values
    lookup_tables['Sigma_eta_13'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_13.csv'), header=None).values
    lookup_tables['Sigma_eta_14'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_14.csv'), header=None).values
    lookup_tables['Sigma_eta_22'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_22.csv'), header=None).values
    lookup_tables['Sigma_eta_23'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_23.csv'), header=None).values
    lookup_tables['Sigma_eta_24'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_24.csv'), header=None).values
    lookup_tables['Sigma_eta_33'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_33.csv'), header=None).values
    lookup_tables['Sigma_eta_34'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_34.csv'), header=None).values
    lookup_tables['Sigma_eta_44'] = pd.read_csv(os.path.join(myDir, 'Sigma_eta_44.csv'), header=None).values
    lookup_tables['mu_X_1'] = pd.read_csv(os.path.join(myDir, 'mu_X_1.csv'), header=None).values
    lookup_tables['mu_X_2'] = pd.read_csv(os.path.join(myDir, 'mu_X_2.csv'), header=None).values
    lookup_tables['Sigma_X_11'] = pd.read_csv(os.path.join(myDir, 'Sigma_X_11.csv'), header=None).values
    lookup_tables['Sigma_X_12'] = pd.read_csv(os.path.join(myDir, 'Sigma_X_12.csv'), header=None).values
    lookup_tables['Sigma_X_22'] = pd.read_csv(os.path.join(myDir, 'Sigma_X_22.csv'), header=None).values
    lookup_tables['mu_CA50_1'] = pd.read_csv(os.path.join(myDir, 'mu_CA50_1.csv'), header=None).values
    lookup_tables['mu_CA50_2'] = pd.read_csv(os.path.join(myDir, 'mu_CA50_2.csv'), header=None).values
    lookup_tables['mu_CA50_3'] = pd.read_csv(os.path.join(myDir, 'mu_CA50_3.csv'), header=None).values
    lookup_tables['mu_CA50_4'] = pd.read_csv(os.path.join(myDir, 'mu_CA50_4.csv'), header=None).values
    lookup_tables['mu_CA50_5'] = pd.read_csv(os.path.join(myDir, 'mu_CA50_5.csv'), header=None).values
    lookup_tables['Sigma_CA50_11'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_11.csv'), header=None).values
    lookup_tables['Sigma_CA50_12'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_12.csv'), header=None).values
    lookup_tables['Sigma_CA50_13'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_13.csv'), header=None).values
    lookup_tables['Sigma_CA50_14'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_14.csv'), header=None).values
    lookup_tables['Sigma_CA50_15'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_15.csv'), header=None).values
    lookup_tables['Sigma_CA50_22'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_22.csv'), header=None).values
    lookup_tables['Sigma_CA50_23'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_23.csv'), header=None).values
    lookup_tables['Sigma_CA50_24'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_24.csv'), header=None).values
    lookup_tables['Sigma_CA50_25'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_25.csv'), header=None).values
    lookup_tables['Sigma_CA50_33'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_33.csv'), header=None).values
    lookup_tables['Sigma_CA50_34'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_34.csv'), header=None).values
    lookup_tables['Sigma_CA50_35'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_35.csv'), header=None).values
    lookup_tables['Sigma_CA50_44'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_44.csv'), header=None).values
    lookup_tables['Sigma_CA50_45'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_45.csv'), header=None).values
    lookup_tables['Sigma_CA50_55'] = pd.read_csv(os.path.join(myDir, 'Sigma_CA50_55.csv'), header=None).values
    lookup_tables['DI_QTY_interp'] = pd.read_csv(os.path.join(myDir, 'DI_QTY_interp.csv'), header=None).values
    lookup_tables['DI_SOI_interp'] = pd.read_csv(os.path.join(myDir, 'DI_SOI_interp.csv'), header=None).values

    return lookup_tables


def load_experiments(myDir, n_cycles_file):

    # Initialize variables
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

    # Path to measured data
    myFiles = [f for f in os.listdir(myDir) if f.endswith('_Data.csv')]

    # Load CSV files
    for fileName in myFiles:

        # Read CSV file
        cycleData = pd.read_csv(os.path.join(myDir, fileName)).values

        # Extract state time series
        M_diesel_exp.append(cycleData[:n_cycles_file, 0])
        M_NH3_exp.append(cycleData[:n_cycles_file, 1])
        M_air_exp.append(cycleData[:n_cycles_file, 2])
        m_diesel_exp.append(cycleData[:n_cycles_file, 3])
        m_NH3_exp.append(cycleData[:n_cycles_file, 4])
        m_air_exp.append(cycleData[:n_cycles_file, 5])
        Q_gross_exp.append(cycleData[:n_cycles_file, 6])
        CA50_exp.append(cycleData[:n_cycles_file, 7])
        DI_SOI_exp.append(cycleData[:n_cycles_file, 8])
        X_res_exp.append(cycleData[:n_cycles_file, 9])

    # Convert to NumPy arrays
    M_diesel_exp = np.concatenate(M_diesel_exp)
    M_NH3_exp = np.concatenate(M_NH3_exp)
    M_air_exp = np.concatenate(M_air_exp)
    m_diesel_exp = np.concatenate(m_diesel_exp)
    m_NH3_exp = np.concatenate(m_NH3_exp)
    m_air_exp = np.concatenate(m_air_exp)
    Q_gross_exp = np.concatenate(Q_gross_exp)
    CA50_exp = np.concatenate(CA50_exp)
    DI_SOI_exp = np.concatenate(DI_SOI_exp)
    X_res_exp = np.concatenate(X_res_exp)
    
    return M_diesel_exp, M_NH3_exp, M_air_exp, \
           m_diesel_exp, m_NH3_exp, m_air_exp, \
           Q_gross_exp, CA50_exp, DI_SOI_exp, X_res_exp


def Gauss_parameters(Diesel_fuel, Diesel_SOI, lookup_tables):
    
    # Extract grid points
    DI_QTY_interp = lookup_tables['DI_QTY_interp'][0, :]
    DI_SOI_interp = lookup_tables['DI_SOI_interp'][:, 0]
    
    # Create the interpolators using RegularGridInterpolator for each lookup table
    interp_mu_eta_1 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_1'].T)
    interp_mu_eta_2 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_2'].T)
    interp_mu_eta_3 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_3'].T)
    interp_mu_eta_4 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_4'].T)
    
    interp_Sigma_eta_11 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_11'].T)
    interp_Sigma_eta_12 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_12'].T)
    interp_Sigma_eta_13 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_13'].T)
    interp_Sigma_eta_14 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_14'].T)
    interp_Sigma_eta_22 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_22'].T)
    interp_Sigma_eta_23 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_23'].T)
    interp_Sigma_eta_24 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_24'].T)
    interp_Sigma_eta_33 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_33'].T)
    interp_Sigma_eta_34 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_34'].T)
    interp_Sigma_eta_44 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_eta_44'].T)
    
    interp_mu_X_1 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_X_1'].T)
    interp_mu_X_2 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_X_2'].T)
    
    interp_Sigma_X_11 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_X_11'].T)
    interp_Sigma_X_12 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_X_12'].T)
    interp_Sigma_X_22 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_X_22'].T)
    
    interp_mu_CA50_1 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50_1'].T)
    interp_mu_CA50_2 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50_2'].T)
    interp_mu_CA50_3 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50_3'].T)
    interp_mu_CA50_4 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50_4'].T)
    interp_mu_CA50_5 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50_5'].T)
    
    interp_Sigma_CA50_11 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_11'].T)
    interp_Sigma_CA50_12 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_12'].T)
    interp_Sigma_CA50_13 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_13'].T)
    interp_Sigma_CA50_14 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_14'].T)
    interp_Sigma_CA50_15 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_15'].T)
    interp_Sigma_CA50_22 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_22'].T)
    interp_Sigma_CA50_23 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_23'].T)
    interp_Sigma_CA50_24 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_24'].T)
    interp_Sigma_CA50_25 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_25'].T)
    interp_Sigma_CA50_33 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_33'].T)
    interp_Sigma_CA50_34 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_34'].T)
    interp_Sigma_CA50_35 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_35'].T)
    interp_Sigma_CA50_44 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_44'].T)
    interp_Sigma_CA50_45 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_45'].T)
    interp_Sigma_CA50_55 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['Sigma_CA50_55'].T)
    
    # Prepare point for interpolation
    point = [Diesel_fuel, Diesel_SOI]

    # Interpolate the values for the given Diesel_fuel and Diesel_SOI
    mu_eta_1_eval = interp_mu_eta_1(point)
    mu_eta_2_eval = interp_mu_eta_2(point)
    mu_eta_3_eval = interp_mu_eta_3(point)
    mu_eta_4_eval = interp_mu_eta_4(point)

    Sigma_eta_11_eval = interp_Sigma_eta_11(point)
    Sigma_eta_12_eval = interp_Sigma_eta_12(point)
    Sigma_eta_13_eval = interp_Sigma_eta_13(point)
    Sigma_eta_14_eval = interp_Sigma_eta_14(point)
    Sigma_eta_22_eval = interp_Sigma_eta_22(point)
    Sigma_eta_23_eval = interp_Sigma_eta_23(point)
    Sigma_eta_24_eval = interp_Sigma_eta_24(point)
    Sigma_eta_33_eval = interp_Sigma_eta_33(point)
    Sigma_eta_34_eval = interp_Sigma_eta_34(point)
    Sigma_eta_44_eval = interp_Sigma_eta_44(point)

    mu_X_1_eval = interp_mu_X_1(point)
    mu_X_2_eval = interp_mu_X_2(point)

    Sigma_X_11_eval = interp_Sigma_X_11(point)
    Sigma_X_12_eval = interp_Sigma_X_12(point)
    Sigma_X_22_eval = interp_Sigma_X_22(point)

    mu_CA50_1_eval = interp_mu_CA50_1(point)
    mu_CA50_2_eval = interp_mu_CA50_2(point)
    mu_CA50_3_eval = interp_mu_CA50_3(point)
    mu_CA50_4_eval = interp_mu_CA50_4(point)
    mu_CA50_5_eval = interp_mu_CA50_5(point)
    
    Sigma_CA50_11_eval = interp_Sigma_CA50_11(point)
    Sigma_CA50_12_eval = interp_Sigma_CA50_12(point)
    Sigma_CA50_13_eval = interp_Sigma_CA50_13(point)
    Sigma_CA50_14_eval = interp_Sigma_CA50_14(point)
    Sigma_CA50_15_eval = interp_Sigma_CA50_15(point)
    Sigma_CA50_22_eval = interp_Sigma_CA50_22(point)
    Sigma_CA50_23_eval = interp_Sigma_CA50_23(point)
    Sigma_CA50_24_eval = interp_Sigma_CA50_24(point)
    Sigma_CA50_25_eval = interp_Sigma_CA50_25(point)
    Sigma_CA50_33_eval = interp_Sigma_CA50_33(point)
    Sigma_CA50_34_eval = interp_Sigma_CA50_34(point)
    Sigma_CA50_35_eval = interp_Sigma_CA50_35(point)
    Sigma_CA50_44_eval = interp_Sigma_CA50_44(point)
    Sigma_CA50_45_eval = interp_Sigma_CA50_45(point)
    Sigma_CA50_55_eval = interp_Sigma_CA50_55(point)

    # Pack outputs
    X_res_mu = np.array([mu_X_1_eval,
                         mu_X_2_eval]).flatten()
    X_res_Sigma = np.array([[Sigma_X_11_eval, Sigma_X_12_eval],
                            [Sigma_X_12_eval, Sigma_X_22_eval]]).reshape(2, 2)

    eta_c_mu = np.array([mu_eta_1_eval,
                         mu_eta_2_eval,
                         mu_eta_3_eval,
                         mu_eta_4_eval]).flatten()
    eta_c_Sigma = np.array([[Sigma_eta_11_eval, Sigma_eta_12_eval, Sigma_eta_13_eval, Sigma_eta_14_eval],
                            [Sigma_eta_12_eval, Sigma_eta_22_eval, Sigma_eta_23_eval, Sigma_eta_24_eval],
                            [Sigma_eta_13_eval, Sigma_eta_23_eval, Sigma_eta_33_eval, Sigma_eta_34_eval],
                            [Sigma_eta_14_eval, Sigma_eta_24_eval, Sigma_eta_34_eval, Sigma_eta_44_eval]]).reshape(4, 4)
    
    CA50_mu = np.array([mu_CA50_1_eval,
                        mu_CA50_2_eval,
                        mu_CA50_3_eval,
                        mu_CA50_4_eval,
                        mu_CA50_5_eval]).flatten()
    CA50_Sigma = np.array([[Sigma_CA50_11_eval, Sigma_CA50_12_eval, Sigma_CA50_13_eval, Sigma_CA50_14_eval, Sigma_CA50_15_eval],
                           [Sigma_CA50_12_eval, Sigma_CA50_22_eval, Sigma_CA50_23_eval, Sigma_CA50_24_eval, Sigma_CA50_25_eval],
                           [Sigma_CA50_13_eval, Sigma_CA50_23_eval, Sigma_CA50_33_eval, Sigma_CA50_34_eval, Sigma_CA50_35_eval],
                           [Sigma_CA50_14_eval, Sigma_CA50_24_eval, Sigma_CA50_34_eval, Sigma_CA50_44_eval, Sigma_CA50_45_eval],
                           [Sigma_CA50_15_eval, Sigma_CA50_25_eval, Sigma_CA50_35_eval, Sigma_CA50_45_eval, Sigma_CA50_55_eval]]).reshape(5, 5)

    return eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, CA50_mu, CA50_Sigma


def conditional_Gauss(mu, Sigma, x2):
    
    # Conditional variance
    sigma_cond = np.sqrt(Sigma[0, 0] - Sigma[0, 1:] @ np.linalg.inv(Sigma[1:, 1:]) @ Sigma[1:, 0])
    # Avoid negative variances
    sigma_cond = max(sigma_cond, 0)

    # Conditional mean
    mu_cond = mu[0] + Sigma[0, 1:] @ np.linalg.inv(Sigma[1:, 1:]) @ (x2 - mu[1:])

    return np.random.normal(mu_cond, sigma_cond)


# Path to lookup tables and experimental data
myDir = 'Dual_fuel_Data/'

# Load lookup tables
lookup_tables = load_lookup_tables(myDir)

# Cycles per file
n_cycles_file = 100

# Load experimental data
M_diesel_exp, M_NH3_exp, M_air_exp, \
m_diesel_exp, m_NH3_exp, m_air_exp, \
Q_gross_exp, CA50_exp, DI_SOI_exp, X_res_exp = load_experiments(myDir, n_cycles_file)

# Simulator constants
AFR_diesel = 14.5 # Stoichiometric air-to-fuel ratio of Diesel (-)
AFR_NH3    = 6.04 # Stoichiometric air-to-fuel ratio of Ammonia (-)
Q_LHV_diesel = 44.1 * 1e6   # Lower heating value of Diesel (J)
Q_LHV_NH3 = 18.6 * 1e6  # Lower heating value of Ammonia (J)
n_cycles = len(Q_gross_exp) # Number of cycles in simulation

# Cost function parameters (NEED TO BE CALIBRATED)
alpha = 8e-06
beta = 0.02
gamma = 0.02

# Scale inputs to simulate how EONS will send messages to ORCAS
# DI_SOI_EONS = (DI_SOI_exp - 42.56) / 6.72
# m_diesel_EONS = (m_diesel_exp * 1e6 - 7.85) / 2.26

# Initialize variables system variables
M_diesel_sim = np.zeros(n_cycles)
M_NH3_sim = np.zeros(n_cycles)
M_air_sim = np.zeros(n_cycles)
eta_c_sim = np.zeros(n_cycles)
Q_gross_sim = np.zeros(n_cycles)
CA50_sim = np.zeros(n_cycles)
X_res_sim = np.zeros(n_cycles)

# Cost function
cost_c = np.zeros(n_cycles)

# Feedback from EONS
Diesel_SOI_EONS = np.zeros(n_cycles)
Diesel_fuel_EONS = np.zeros(n_cycles)

# Initial condition
M_diesel_sim[0] = M_diesel_exp[0]
M_NH3_sim[0] = M_NH3_exp[0]
M_air_sim[0] = M_air_exp[0]

lll = 1000
# Loop through cycles
# for i in range(n_cycles-1):

for i in range(lll):

    # Diesel feedback inputs: THESE ARE THE INPUTS YOU CAN MODIFY
    Diesel_SOI_EONS[i]  = 0
    Diesel_fuel_EONS[i] = 0

    # Diesel, Ammonia, and Air feedforward inputs:
    # THESE WILL BE PROVIDED FROM ORCAS
    Diesel_SOI_FF = DI_SOI_exp[i]             # deg bTDC
    Diesel_fuel_FF = m_diesel_exp[i] * 1e6    # mg
    Ammonia_fuel_FF = m_NH3_exp[i] * 1e6      # mg
    Fresh_air_FF = m_air_exp[i] * 1e6         # mg

    # Combine EONS feedback and feelforward diesel commands
    Diesel_SOI_total =  Diesel_SOI_FF + Diesel_SOI_EONS[i]
    Diesel_fuel_total = Diesel_fuel_FF + Diesel_fuel_EONS[i]

    # Gaussian parameters
    eta_c_mu, eta_c_Sigma, X_res_mu, X_res_Sigma, CA50_mu, CA50_Sigma \
        = Gauss_parameters(Diesel_fuel_total, Diesel_SOI_total, lookup_tables)
    
    # Mass state
    state = np.array([M_diesel_sim[i], M_NH3_sim[i], M_air_sim[i]])
    # Combustion state
    state_CA50 = np.array([M_diesel_sim[i] * 1e6,
                           M_NH3_sim[i] * 1e6,
                           M_air_sim[i] * 1e4,
                           Diesel_SOI_total])
    
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
                                [      1 - eta_c_sim[i],                0,            0],
                                [             0,                1 -  eta_c_sim[i],    0],
                                [-AFR_diesel * eta_c_sim[i], -AFR_NH3 * eta_c_sim[i], 1]])

    # Fresh fuel and air (Kg)
    input_mass = np.array([Diesel_fuel_total,
                           Ammonia_fuel_FF,
                           Fresh_air_FF]) * 1e-6

    # Dynamics
    if i < n_cycles - 1:
        next_state = Matrix_res @ state + input_mass
        M_diesel_sim[i+1] = next_state[0]
        M_NH3_sim[i+1] = next_state[1]
        M_air_sim[i+1] = next_state[2]

#%% Plot 
plt.figure()
plt.plot(M_diesel_exp, '-o', label='Experimental')
plt.plot(M_diesel_sim, '-o', label='Simulated')
plt.ylabel('In-Cylinder diesel (kg)')
plt.legend()
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

plt.figure()
plt.plot(M_NH3_exp, '-o', label='Experimental')
plt.plot(M_NH3_sim, '-o', label='Simulated')
plt.ylabel('In-Cylinder diesel (kg)')
plt.legend()
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

plt.figure()
plt.plot(M_air_exp, '-o', label='Experimental')
plt.plot(M_air_sim, '-o', label='Simulated')
plt.ylabel('In-Cylinder diesel (kg)')
plt.legend()
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

# Plot Q_gross
plt.figure()
plt.plot(Q_gross_exp, '-o', label='Experimental')
plt.plot(Q_gross_sim, '-o', label='Simulated')
plt.ylabel('Q_gross (J)')
plt.legend()
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

plt.figure()
plt.plot(eta_c_sim, '-o', label='Simulated')
plt.legend()
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

# Plot CA50
plt.figure()
plt.plot(CA50_exp, label='Experimental')
plt.plot(CA50_sim, label='Simulated')
plt.ylabel('CA50 (aTDC)')
plt.legend()
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

# Plot combustion efficiency
plt.figure()
plt.plot(eta_c_sim)
plt.legend()
plt.title('Combustion efficiency')
plt.xlabel('Cycles')
plt.ylim([0, 1])
plt.xlim([95, lll])
plt.show()

# Plot cost function
plt.figure()
plt.plot(cost_c)
plt.legend()
plt.title('Cost function')
plt.xlabel('Cycles')
plt.xlim([95, lll])
plt.show()

# Plot feedback from EONS
fig, ax = plt.subplots(2, 1)
ax[0].plot(Diesel_SOI_EONS)
ax[0].set_title('Feedback from EONS')
ax[0].set_ylabel('Diesel SOI')
ax[1].plot(Diesel_fuel_EONS)
plt.xlim([95, lll])
ax[1].set_ylabel('Diesel mass')
ax[1].set_xlabel('Cycle')

# Plot feedback from EONS
fig, ax = plt.subplots(3, 1)
ax[0].plot(M_diesel_sim)
plt.xlim([95, lll])
ax[1].plot(M_NH3_sim)
ax[2].plot(M_air_sim)