import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

def load_lookup_tables(myDir):
    
    # Load lookup tables from files
    lookup_tables = {}
    lookup_tables['mu_eta_1'] = pd.read_csv(os.path.join(myDir, 'mu_eta_1.csv'), header=None).values
    lookup_tables['mu_CA50'] = pd.read_csv(os.path.join(myDir, 'mu_CA50.csv'), header=None).values
    lookup_tables['DI_QTY_interp'] = pd.read_csv(os.path.join(myDir, 'DI_QTY_interp.csv'), header=None).values
    lookup_tables['DI_SOI_interp'] = pd.read_csv(os.path.join(myDir, 'DI_SOI_interp.csv'), header=None).values

    return lookup_tables


def average_cost(point, alpha, lookup_tables):

    # Unpack variables
    Diesel_fuel, Diesel_SOI = point
    
    # Extract grid points
    DI_QTY_interp = lookup_tables['DI_QTY_interp'][0, :]
    DI_SOI_interp = lookup_tables['DI_SOI_interp'][:, 0]
    
    # Average combustion efficiency
    interp_mu_eta = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_eta_1'].T)
    
    # Average CA50
    interp_mu_CA50 = RegularGridInterpolator((DI_QTY_interp, DI_SOI_interp), lookup_tables['mu_CA50'].T)

    # Interpolate the values for the given Diesel_fuel and Diesel_SOI
    mu_eta  = interp_mu_eta(point) / 100
    mu_CA50 = interp_mu_CA50(point)
    
    # Cost function
    cost_c = alpha * (mu_CA50 - 7.5) ** 2 + (mu_eta - 1) ** 2

    return cost_c


# Path to measured data
myDir = 'Model_Data/'
myFiles = [f for f in os.listdir(myDir) if f.endswith('_Data.csv')]

# Load lookup tables
lookup_tables = load_lookup_tables(myDir)

# Initial guess
point_0 = [np.mean(lookup_tables['DI_QTY_interp'][0, :]), np.mean(lookup_tables['DI_SOI_interp'][:, 0])]

# Bounds
bounds = ((np.min(lookup_tables['DI_QTY_interp'][0, :]), np.max(lookup_tables['DI_QTY_interp'][0, :])),
          (np.min(lookup_tables['DI_SOI_interp'][:, 0]), np.max(lookup_tables['DI_SOI_interp'][:, 0])))

# Cost function parameter
alpha = 0.02

# Optimize
point_opt = minimize(average_cost,
                     point_0,
                     args=(alpha, lookup_tables),
                     method='SLSQP',
                     jac='3-point',
                     bounds=bounds,
                     options={'disp': False}).x

# Cost
cost_opt = average_cost(point_opt, alpha, lookup_tables)

# Solution
print(f"Point = [{point_opt[0]:.2f}, {point_opt[1]:.2f}], cost = {cost_opt[0]:.2f}")