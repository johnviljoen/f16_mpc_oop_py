#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 22:02:06 2020

@author: johnviljoen
"""

''' This file contains all parameters (paras) required to run the simulation,
the aircraft, environmental, simulation, initial conditions, and other parameters'''

#import numpy as np
import numpy as np
from numpy import pi
from scipy.constants import g

# In[simulation parameters]  

time_step, time_start, time_end = 0.001, 0., 3.

# fi_flag = 1 -> high fidelity model (full Nguyen)
# fi_flag = 1 -> low fidelity model (Stevens Lewis reduced)
fi_flag = 1

# stability_flag only functional for high fidelity model currently!
# stability_flag = 1 -> unstable xcg 35% model
# stability_flag = 0 -> stable xcg 25% model
stab_flag = 0

# In[MPC parameters]

hzn = 10

pred_dt = 0.001

# In[initial_conditions]  

''' states in m/s, rad, rad/s '''
npos        = 0.                # m
epos        = 0.                # m
h           = 3048.             # m
phi         = 0.                # rad
theta       = 0.                # rad
psi         = 0.                # rad

vt          = 213.36            # m/s
alpha       = 1.0721 * pi/180   # rad
beta        = 0.                # rad
p           = 0.                # rad/s
q           = 0.                # rad/s
r           = 0.                # rad/s

''' control states in lbs, deg '''
T           = 2886.6468         # lbs
dh          = -2.0385           # deg
da          = -0.087577         # deg
dr          = -0.03877          # deg
lef         = 0.3986            # deg

# In[limits]

npos_min        = -np.inf       # (m)
epos_min        = -np.inf       # (m)
h_min           = 0             # (m)
phi_min         = -np.inf       # (deg)
theta_min       = -np.inf       # (deg)
psi_min         = -np.inf       # (deg)
V_min           = 0             # (m/s)
alpha_min       = -20.          # (deg)
beta_min        = -30.          # (deg)
p_min           = -30           # (deg/s)
q_min           = -10           # (deg/s)
r_min           = -5            # (deg/s)

T_min           = 1000          # (lbs)
dh_min          = -25           # (deg)
da_min          = -21.5         # (deg)
dr_min          = -30.          # (deg)
lef_min         = 0.            # (deg)

npos_max        = np.inf        # (m)
epos_max        = np.inf        # (m)
h_max           = 10000         # (m)
phi_max         = np.inf        # (deg)
theta_max       = np.inf        # (deg)
psi_max         = np.inf        # (deg)
V_max           = 900           # (m/s)
alpha_max       = 90            # (deg)
beta_max        = 30            # (deg)
p_max           = 30            # (deg/s)
q_max           = 10            # (deg/s)
r_max           = 5             # (deg/s)

T_max           = 19000         # (lbs)
dh_max          = 25            # (deg)
da_max          = 21.5          # (deg)
dr_max          = 30            # (deg)
lef_max         = 25            # (deg)

# In[wrap for input]  

# initial_state_vector = np.array([npos, epos, h, phi, theta, psi, vt, alpha, beta, p, q, r, T, dh, da, dr, lef, fi_flag])

simulation_parameters = [time_step, time_start, time_end, stab_flag, fi_flag]

paras_mpc = [hzn, pred_dt]

m2f = 3.28084 # metres to feet conversion
f2m = 1/m2f # feet to metres conversion
initial_state_vector_ft_rad = np.array([npos*m2f, epos*m2f, h*m2f, phi, theta, psi, vt*m2f, alpha, beta, p, q, r, T, dh, da, dr, lef, -alpha*180/pi])
    
act_lim = [[T_max, dh_max, da_max, dr_max, lef_max],
           [T_min, dh_min, da_min, dr_min, lef_min]]

x_lim = [[npos_max, epos_max, h_max, phi_max, theta_max, psi_max, V_max, alpha_max, beta_max, p_max, q_max, r_max],
         [npos_min, epos_min, h_min, phi_min, theta_min, psi_min, V_min, alpha_min, beta_min, p_min, q_min, r_min]]

# In[additional info provided for brevity]

# weight                  = 91188         # Newtons

# Ixx                     = 12875         # Kg m^2
# Iyy                     = 75674         # Kg m^2
# Izz                     = 85552         # Kg m^2
# Ixz                     = 1331          # Kg m^2
# # the other Izy, Iyz = 0

# b                       = 9.144         # m wingspan
# S                       = 27.87         # m^2 wing area
# cbar                    = 3.45          # m wing mean aerodynamic chord

# He                      = 216.9         # engine angular momentum constant

# x_cg_ref                = 0.35 * cbar   # assuming mac = cbar
# x_cg                    = 0.8*x_cg_ref  # FOR NOW THIS IS WRONG

# # unecessary:
# length = 14.8 #m
# height = 4.8 #m